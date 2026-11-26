%% Setup: Adaptive MPC + Recursive ARX for PEMFC + Battery SoC (target 0.92)
clear; clc;

%% 1. Dimensions and timing
Ts = 60; % sec
nu  = 5; % MVs: [current_density, pressure, lambda, temperature, humidity]
nmd = 1; % MD : [load]
ny  = 1;% OV : [SoC]

%% 2. Nominal operating point (physical units, inside the MV bounds)
% order: [current_density, pressure, lambda, temperature, humidity, load]
u0 = [0.80; 150; 1.7; 70; 0.35; 0]; % MVs + MD
y0 = 0.92; % SoC operating point (= target)
x0 = 0.92; % single plant state = SoC (C0 = 1)

%% 3. Seed plant model (discrete SS)
% Physics: SoC integrates net charging power, so the seed is an integrator
% driven mainly by current_density (+) and load (-); ARX learns the rest.
A0 = 1; % SoC integrator
% j, p, lam, T, RH, load
B0 = [0.0020, 0.0003, 0.0002, 0.0003, 0.0002, -0.0010] * Ts;
C0 = 1;
D0 = zeros(ny, nu + nmd);

plant0 = ss(A0, B0, C0, D0, Ts);
plant0.InputName  = {'current_density','pressure','lambda','temperature','humidity','load'};
plant0.OutputName = {'SoC'};
plant0 = setmpcsignals(plant0, 'MV', 1:nu, 'MD', nu + (1:nmd));

%% 4. MPC controller
p = 20;  m = 4; % prediction / control horizons
mpcobj = mpc(plant0, Ts, p, m);

mpcobj.Model.Nominal.U = u0; % feasible nominal (inside MV bounds)
mpcobj.Model.Nominal.Y = y0;
mpcobj.Model.Nominal.X = x0;
mpcobj.Model.Nominal.DX = 0; % one state

% Output: keep SoC near 0.92
mpcobj.OV(1).Min = 0.80; 
mpcobj.OV(1).Max = 1.00;
mpcobj.Weights.OutputVariables = 10;

% MV bounds (physical units), move-rate limits, per channel
mv_min = [0.30,  50, 1.6, 55, 0.30];
mv_max = [1.80, 200, 2.1, 75, 0.90];
mv_rate = [0.10,  25, 0.25, 2.0, 0.05];
for k = 1:nu
    mpcobj.MV(k).Min = mv_min(k);
    mpcobj.MV(k).Max = mv_max(k);
    mpcobj.MV(k).RateMin = -mv_rate(k);
    mpcobj.MV(k).RateMax = mv_rate(k);
end

% Weights: current_density drives SoC (free); anchor the weak MVs near nominal
mpcobj.Weights.ManipulatedVariables = [0, 0.05, 0.05, 0.05, 0.05];
mpcobj.Weights.ManipulatedVariablesRate = 0.1 * ones(1, nu);

% Scale factors = each signal's span (MVs/MD/OV differ by orders of magnitude)
sf_mv = [1.5, 150, 0.5, 20, 0.6];
for k = 1:nu, mpcobj.MV(k).ScaleFactor = sf_mv(k); end
mpcobj.DV(1).ScaleFactor = 1; % set to your load range
mpcobj.OV(1).ScaleFactor = 0.2; % SoC operating span (0.8..1.0)

%% 5. MPC initial controller state (mpcstate)
soc0 = 0.30; % actual battery SoC at t = 0
xmpc0 = mpcstate(mpcobj);
xmpc0.Plant = soc0; % C0=1 -> state == SoC
xmpc0.LastMove = u0(1:nu)'; % last MV, length nu

%% 6. Recursive ARX estimator
% seed via InitialA / InitialB
na = 1;
nb = ones(1, nu + nmd); % one b-coeff per input (5 MV + 1 MD)
nk = ones(1, nu + nmd); % one-step input delay per input
rARX = recursiveARX([na nb nk]);
rARX.ForgettingFactor = 0.99;

% recursiveARX has NO InitialParameters field. Seed A(q) and B(q) instead:
rARX.InitialA = [1 -1]; % length na+1=3, first elem must be 1
b_seed = [0.0020; 0.0003; 0.0002; 0.0003; 0.0002; -0.0010];   % j, p, lam, T, RH, load
rARX.InitialB = [zeros(nu+nmd,1), b_seed]; % Nu-by-max(nb+nk) = 6-by-2, each row [0  b_i]
rARX.InitialParameterCovariance = 1e2; % < default 1e4 so the seed carries some weight