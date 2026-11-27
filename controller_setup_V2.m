%% Note on improvements
%{
adding temperature control to make use of some other controller devices
and adding another target of minimizing the voltage loss (mitigate
degradation). Guessing the prediction should be able to choose the option
that slow down the aging most, according to the degradation function
(aka. improvised penalty function - voltage loss) while regulating both
SoC and cell temperature.
%}

%% Setup: Adaptive MPC + Recursive ARX for PEMFC — SoC, Cell Temp, Degradation (MIMO)
% OVs: [SoC (target 0.92), cell temperature (target 70C), voltage-loss rate (target 0)]
% MVs: [current_density, pressure, lambda, cooling, humidity]   MD: [load]
% Only apply for simulate_model.slx
clear all; clc;

%% 1. Dimensions and timing
Ts  = 60;   % sec
nu  = 5;    % MVs: [current_density, pressure, lambda, cooling, humidity]
nmd = 1;    % MD : [load]
ny  = 3;    % OVs: [SoC, cell temperature, voltage-loss rate]

%% 2. Nominal operating point (physical units, inside the MV bounds)
% order: [current_density, pressure, lambda, cooling, humidity, load]
u0 = [0.80; 150; 1.7; 0.50; 0.35; 0];   % MVs + MD
y0 = [0.92; 70; 0];                     % [SoC; T_cell; Vloss_rate] operating point
x0 = [0.92; 70; 0];                     % states = outputs (C0 = I)

%% 3. Seed plant model (discrete SS, 3 states)
% SoC    : integrator driven by current_density (+) and load (-).
% T_cell : first-order thermal, heated by current_density & load, cooled by 'cooling'.
% Vloss  : instantaneous voltage-loss (degradation) rate ~ static map of operation.
A0 = [1   0    0;        % SoC integrator
      0   0.9  0;        % T_cell first-order pole (tau ~ Ts/(1-0.9) = 600 s)
      0   0    0];       % Vloss_rate: static (no memory)

% columns: j, p, lam, cooling, RH, load
B0 = zeros(ny, nu+nmd);
B0(1,:) = [0.0020, 0.0003, 0.0002, 0.0000, 0.0002, -0.0010] * Ts;  % SoC   (per-step)
B0(2,:) = [1.00,   0.00,  -0.05,  -1.00,   0.00,    0.50];          % T_cell (per-step)
B0(3,:) = [0.05,   0.00,   0.00,  -0.01,  -0.02,    0.01];          % Vloss rate
%          ^ rises with current density & load, falls with cooling & humidity.
%          Swap these for your fitted cathode-degradation regression gains.

C0 = eye(ny);
D0 = zeros(ny, nu + nmd);

plant0 = ss(A0, B0, C0, D0, Ts);
plant0.InputName  = {'current_density','pressure','lambda','cooling','humidity','load'};
plant0.OutputName = {'SoC','T_cell','Vloss_rate'};
plant0 = setmpcsignals(plant0, 'MV', 1:nu, 'MD', nu + (1:nmd));

%% 4. MPC controller
p = 20;  m = 4;                 % prediction / control horizons
w_deg = 1;                      % degradation vs performance tradeoff knob (raise = gentler)
mpcobj = mpc(plant0, Ts, p, m);
mpcobj.Model.Nominal.U  = u0;   % feasible nominal (inside MV bounds)
mpcobj.Model.Nominal.Y  = y0;
mpcobj.Model.Nominal.X  = x0;
mpcobj.Model.Nominal.DX = [0; 0; 0];

% Output bounds (soft) and tracking priority
mpcobj.OV(1).Min = 0.80;  mpcobj.OV(1).Max = 1.00;   % SoC
mpcobj.OV(2).Min = 60;    mpcobj.OV(2).Max = 85;     % T_cell (C)
mpcobj.OV(3).Min = 0;     mpcobj.OV(3).Max = Inf;    % Vloss rate (>= 0)
mpcobj.Weights.OutputVariables = [10, 5, w_deg];     % [SoC, T_cell, Vloss]

% MV bounds (physical units), move-rate limits, per channel
%         j     p    lam  cool  RH
mv_min  = [0.30,  50, 1.6, 0.10, 0.30];
mv_max  = [1.80, 200, 2.1, 1.00, 0.90];
mv_rate = [0.10,  25, 0.25,0.10, 0.05];
for k = 1:nu
    mpcobj.MV(k).Min     =  mv_min(k);
    mpcobj.MV(k).Max     =  mv_max(k);
    mpcobj.MV(k).RateMin = -mv_rate(k);
    mpcobj.MV(k).RateMax =  mv_rate(k);
end

% Weights: current_density free (drives SoC), cooling free (drives T_cell),
% anchor the weak MVs (p, lam, RH) near nominal.
%                                    j   p    lam  cool  RH
mpcobj.Weights.ManipulatedVariables     = [0, 0.05, 0.05, 0,  0.05];
mpcobj.Weights.ManipulatedVariablesRate = 0.1 * ones(1, nu);

% Scale factors = each signal's span
sf_mv = [1.5, 150, 0.5, 1.0, 0.6];      % j, p, lam, cooling, RH
for k = 1:nu, mpcobj.MV(k).ScaleFactor = sf_mv(k); end
mpcobj.DV(1).ScaleFactor = 1;           % set to your load range
mpcobj.OV(1).ScaleFactor = 0.2;         % SoC span (0.8..1.0)
mpcobj.OV(2).ScaleFactor = 20;          % T_cell span (~60..85 C)
mpcobj.OV(3).ScaleFactor = 0.1;         % Vloss-rate span (tune to your units)

%% 5. MPC initial controller state (mpcstate)
soc0 = 0.30;   % actual battery SoC at t = 0
T0   = 40;     % actual cell temperature at t = 0 (cold start, C)
V0   = 0;      % voltage-loss rate at t = 0
xmpc0 = mpcstate(mpcobj);
xmpc0.Plant    = [soc0; T0; V0];   % C0 = I -> state == output
xmpc0.LastMove = u0(1:nu)';        % last MV, length nu

%% 6. Recursive ARX estimator — bank of first-order MISO estimators (one per output)
% Each output is MISO over the 6 inputs, seeded like the original SISO estimator.
% Stacking gives the MIMO prediction model (diagonal A since na=1, nk=1):
%     x(k+1) = diag(-a_i) x(k) + [b_1; b_2; b_3] u(k),   y = x
na = 1;
nb = ones(1, nu + nmd);     % one b-coeff per input
nk = ones(1, nu + nmd);     % one-step input delay per input

% A(q) seeds:  SoC = integrator,  T_cell = stable pole 0.9,  Vloss = static
a_seed = { [1 -1], [1 -0.9], [1 0] };
% B(q) seeds (per input, no Ts):            j        p        lam     cool     RH      load
b_seed = { [0.0020; 0.0003; 0.0002; 0.0000; 0.0002; -0.0010];   % SoC
           [1.00;   0.00;  -0.05;  -1.00;   0.00;    0.50   ];   % T_cell
           [0.05;   0.00;   0.00;  -0.01;  -0.02;    0.01   ] }; % Vloss rate

rARX = cell(ny,1);
for i = 1:ny
    e = recursiveARX([na nb nk]);
    e.ForgettingFactor            = 0.99;
    e.InitialA                    = a_seed{i};                    % first elem must be 1
    e.InitialB                    = [zeros(nu+nmd,1), b_seed{i}];  % 6x2, each row [0 b_i]
    e.InitialParameterCovariance  = 1e2;    % < default 1e4 so the seed carries weight
    rARX{i} = e;
end

%% --- Reference & assembling the MIMO model for the Adaptive MPC block ---
% Reference (Simulink ref port or mpcmoveAdaptive):  ref = [0.92, 70, 0]
% One Recursive Polynomial Model Estimator block per output, then a MATLAB Function:
%     A = diag([-a1(2); -a2(2); -a3(2)]);   % -a_i coefficient per output
%     B = [b1; b2; b3];                      % each b_i is 1 x (nu+nmd)
%     C = eye(ny);  D = zeros(ny, nu+nmd);
%     sys = setmpcsignals(ss(A,B,C,D,Ts), 'MV',1:nu, 'MD',nu+(1:nmd));
%  -> feed sys.A, sys.B (+ nominal) to the Adaptive MPC block's model port.
