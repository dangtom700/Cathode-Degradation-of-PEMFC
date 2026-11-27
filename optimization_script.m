clear; clc;
age_test = [0, 1150, 2300];
% order: [current_density, pressure, lambda, T_in, humidity]
lb = [0.3  50 1.6 55 0.30];
ub = [1.8 200 2.1 75 0.90];
w = 0.5; % Weight for voltage loss

rng(0,"twister");
num_samples = 10;
r = lb + (ub-lb).*rand(num_samples, numel(lb));   % 10x5 independent starts

o1 = optimoptions('fmincon','Algorithm','sqp','Display','off');
o2 = optimoptions('ga','Display','off','PopulationSize',60);

rows = [];
for age_h = age_test
    for k = 1:num_samples
        x0 = r(k,:);
        [xf,Jf,ef] = fmincon(@(x)neg_netpower(x,age_h,w),x0,[],[],[],[],lb,ub,...
                             @(x)netpower_con(x,age_h),o1);
        [xg,Jg,eg] = ga(@(x)neg_netpower(x,age_h,w),numel(lb),[],[],[],[],lb,ub,...
                        @(x)netpower_con(x,age_h),o2);
        rows = [rows; age_h, k, -Jf, xf, ef, -Jg, xg, eg];
    end
end

vn = {'age_h','start','fmincon_Pnet','f_j','f_pc','f_lam','f_Tin','f_RH','f_flag',...
      'GA_Pnet','g_j','g_pc','g_lam','g_Tin','g_RH','g_flag'};
T = array2table(rows,'VariableNames',vn);
T = varfun(@(c)round(c,3), T);  T.Properties.VariableNames = vn;   % tidy
disp(T)
writetable(T,'optim_runs.csv');   % for the report

%% ---- objective & constraints ----
function J = neg_netpower(x, age_h, w)

    j=x(1); p_c=x(2); lam=x(3); T_in=x(4); RH=x(5);
    T_cell = cell_temperature(j, age_h, T_in);
    Pstack = cellmodel(j, lam, age_h, T_cell, p_c, T_in, RH);

    % Normalized
    loss = -voltage_loss(j, p_c, lam, T_in, RH, age_h);
    u_vector = [Pstack loss] / norm([Pstack loss]);
    J = -[1-w w] * u_vector';
end

function [c, ceq] = netpower_con(x, age_h)
    j=x(1); p_c=x(2); lam=x(3); T_in=x(4); RH=x(5);
    T_cell = cell_temperature(j, age_h, T_in);
    Pstack = cellmodel(j, lam, age_h, T_cell, p_c, T_in, RH);
    U = Pstack / (j*280*10); % back out per-cell voltage from power

    c   = [ T_cell - 80; % cell temperature <= 80 C
            0.6 - U]; % U >= 0.6 V
    ceq = [];
end

%% Cell temperature and voltage
function T = cell_temperature(current_density, age_hour, T_in)
    arguments
        current_density double {mustBeNonnegative}
        age_hour double {mustBeNonnegative}
        T_in double = 55
    end
    
    a = -0.3687; b = 1.6828; c = 0.2772; d = -0.0154; e = -0.1059;
    F_cl = 6.0; % coolant flow [L/min] (fixed)
    
    T = T_in + a + b.*current_density.*6./F_cl ...
        + c.*current_density.^2 + d.*(T_in - 70) ...
        + e.*age_hour./1000;
end

function cell_power = cellmodel(j, lam, age_h, T_cell, p_c, T_in, RH)
    arguments
        j double {mustBeNonnegative}
        lam double {mustBePositive}
        age_h double {mustBeNonnegative}
        T_cell double
        p_c double = 150
        T_in double = 55
        RH double {mustBeGreaterThanOrEqual(RH,0),mustBeLessThanOrEqual(RH,1)} = 0.35
    end
    psat = @(t) 0.61094 .* exp(17.625 .* t ./ (t + 243.04));
    
    % physical constants
    R = 8.314; F = 96485; P0 = 101.325; % J/mol/K, C/mol, kPa
    Ea = 80; jx = exp(-6); tau = 300; % kJ/mol, A/cm^2, h
    % fitted voltage-model coefficients
    lnj0 = -13.6563; j0f = -0.4934; j0s = 0.0801;
    alpha = 0.9740;  gam = 0.4081;  nu  = 0.1985;
    lnc  = -6.9859;  m   = 2.0241;  cf  = 0.0897; cs = 0.2336; kT = 0.0117;
    hb0  = 0.0483;   hb1 = 0.0018;  hb2 = 0.0038; % HFR = hb0 + hb1*age/1000 + hb2*j
    % secondary operating params (fixed)
    lam_a = 1.5; RH_a = 0.50; dpa = 10; % anode stoich, anode RH, anode-cathode p offset
    cell_area = 280; num_cell = 10;
    
    T = T_cell + 273.15; % K
    lam = max(lam, 1.001);
    pa  = p_c + dpa; % anode gauge pressure [kPaG]
    
    % water vapour partial pressures (RH referenced to coolant inlet)
    pw_c = RH   .* psat(T_in);
    pw_a = RH_a .* psat(T_in);
    
    % reactant partial pressures (log-mean O2 along the cathode channel) [kPa]
    pO2 = 0.2095 .* (p_c + P0 - pw_c) .* (1./lam) ./ (-log(1 - 1./lam));
    pH2 = 0.80   .* (pa  + P0 - pw_a) .* (1 - 1./(2.*lam_a));
    
    % age factors
    gf = 1 - exp(-age_h ./ tau);   gs = age_h ./ 1000;
    
    % Nernst potential
    E = 1.229 - 8.46e-4 .* (T - 298.15) + R.*T./(2.*F) .* log(pH2./P0 .* sqrt(pO2./P0));
    
    % activation (ORR): exchange current drops with age
    j0  = exp(lnj0 + j0f.*gf + j0s.*gs) .* exp(-Ea.*1e3./R .* (1./T - 1./343.15));
    act = R.*T./(alpha.*F) .* (log((j + jx)./j0) - gam.*log(pO2./P0));
    
    % ohmic / membrane: HFR rises linearly with age (imposed from EIS)
    HFR = hb0 + hb1.*age_h./1000 + hb2.*j;
    ohm = j .* HFR .* (RH./0.35).^(-nu);
    
    % mass transport: amplitude grows with age
    mt = exp(lnc) .* (1 + cf.*gf + cs.*gs) .* exp(-kT.*(T - 343.15)) .* (j ./ (pO2./P0)).^m;
    
    % Final voltage
    U = E - act - ohm - mt;
    
    % Power
    cell_power = U .* j .* cell_area .* num_cell; % W
end

function loss = voltage_loss(j, p_c, lam, T_in, RH, age_h)
    % Per-cell aging voltage loss vs BoL at the SAME setpoint [mV]
    A=280; N=10;
    T_age = cell_temperature(j, age_h, T_in);
    T_bol = cell_temperature(j, 0,     T_in); % BoL thermal at same setpoint
    U_age = cellmodel(j, lam, age_h, T_age, p_c, T_in, RH)/(j*A*N);
    U_bol = cellmodel(j, lam, 0,     T_bol, p_c, T_in, RH)/(j*A*N);
    loss = U_bol - U_age; % V per cell
end

