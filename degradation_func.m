function out = degradation_func(j, age_h, opts)
%DEGRADATION_FUNC  Semi-empirical PEMFC cell-voltage / degradation model.
%   Single-call version of the voltage model fitted in
%   optimize_operating_point.m (local function vparts). Given a load and a
%   stack age it returns the cell voltage, the Nernst potential, the three
%   loss terms, the mean cell temperature, and the voltage lost to aging.
%
%   Degradation enters through two age factors (t = age in hours):
%       g_f = 1 - exp(-t/tau)      fast break-in   (tau = 300 h)
%       g_s = t/1000               slow linear wear
%   acting on the ORR exchange current (j0), the membrane resistance (HFR)
%   and the mass-transport amplitude.
%
%   USAGE
%       out = degradation_func(j, age_h)                 % nominal setpoint
%       out = degradation_func(j, age_h, Name=Value)     % custom setpoint
%
%   INPUTS
%       j       current density [A/cm^2]   (scalar or array)
%       age_h   stack age [h]              (scalar, or array matching j)
%
%   NAME-VALUE OPTIONS (defaults = nominal operating point)
%       p_c      cathode OUTLET pressure [kPaG]          (150)
%       lambda_c cathode stoichiometry [-]               (1.7)
%       T_in     coolant inlet temperature [degC]        (70)
%       RH_c     cathode inlet RH, ref. coolant inlet    (0.35)
%       F_cl     coolant flow [L/min]                    (6.0)
%       lambda_a anode stoichiometry [-]                 (1.5)
%       RH_a     anode inlet RH, ref. coolant inlet      (0.50)
%       dpa      anode-cathode pressure offset [kPa]     (10)
%       T_cell   mean cell temperature [degC] (override; [] = use model)
%       params   struct of fitted coefficients (override the defaults)
%
%   OUTPUT (struct, element-wise over j / age_h)
%       .U            cell voltage [V]
%       .E            Nernst potential [V]
%       .eta_act      activation (ORR) loss [V]
%       .eta_ohm      ohmic / membrane loss [V]
%       .eta_mt       mass-transport loss [V]
%       .T_cell       mean cell temperature [degC]
%       .j0           ORR exchange current density [A/cm^2]
%       .HFR          high-frequency resistance [Ohm*cm^2]
%       .dU_aging_mV  voltage lost vs BoL at the same setpoint [mV]
%       .params       the coefficient struct actually used
%
%   NOTE  The HFR age slope is imposed from Tier-1 EIS, not fitted here.
%         Trust the condition slopes (gamma, nu, kT, m) only inside the
%         tested range: p_c 50-200 kPaG, lambda_c 1.6-2.1, T 55-75 degC,
%         RH 0.30-0.90, j up to ~2.3 A/cm^2.

    arguments
        j        double
        age_h    double
        opts.p_c      double = 150
        opts.lambda_c double = 1.7
        opts.T_in     double = 70
        opts.RH_c     double = 0.35
        opts.F_cl     double = 6.0
        opts.lambda_a double = 1.5
        opts.RH_a     double = 0.50
        opts.dpa      double = 10
        opts.T_cell   double = []
        opts.params   struct = default_params()
    end

    P = opts.params;

    % --- mean cell temperature (semi-empirical thermal model) ------------
    if isempty(opts.T_cell)
        t = P.th_T;    % [a, b_flow, c_j2, d_Tin, e_age]
        T_cell = opts.T_in + t(1) + t(2).*j.*6./opts.F_cl + t(3).*j.^2 ...
                 + t(4).*(opts.T_in - 70) + t(5).*age_h./1000;
    else
        T_cell = opts.T_cell;
    end

    % --- voltage and loss breakdown at this age -------------------------
    [U, E, act, ohm, mt, j0, HFR] = vcell(P, j, age_h, T_cell, opts);

    % --- voltage lost to aging: same setpoint, age = 0 ------------------
    if isempty(opts.T_cell)
        T_bol = opts.T_in + P.th_T(1) + P.th_T(2).*j.*6./opts.F_cl ...
                + P.th_T(3).*j.^2 + P.th_T(4).*(opts.T_in - 70);
    else
        T_bol = opts.T_cell;
    end
    U_bol = vcell(P, j, 0, T_bol, opts);

    out = struct('U', U, 'E', E, 'eta_act', act, 'eta_ohm', ohm, ...
                 'eta_mt', mt, 'T_cell', T_cell, 'j0', j0, 'HFR', HFR, ...
                 'dU_aging_mV', 1e3.*(U_bol - U), 'params', P);
end

function [U, E, act, ohm, mt, j0, HFR] = vcell(P, j, age_h, T_cell, opts)
    R = P.R; F = P.F; P0 = P.P0;
    T = T_cell + 273.15;                      % K

    % dew points from RH referenced to coolant inlet (as in fcSystem)
    dptc = dewpoint(opts.RH_c .* psat(opts.T_in));
    dpta = dewpoint(opts.RH_a .* psat(opts.T_in));
    pc = opts.p_c;  pa = opts.p_c + opts.dpa; % kPaG
    lam  = max(opts.lambda_c, 1.001);
    lama = opts.lambda_a;

    % partial pressures (log-mean O2 along the cathode channel), kPa
    pO2 = 0.2095 .* (pc + P0 - psat(dptc)) .* (1./lam) ./ (-log(1 - 1./lam));
    pH2 = 0.80  .* (pa + P0 - psat(dpta)) .* (1 - 1./(2.*lama));

    % age factors
    gf = 1 - exp(-age_h./P.tau);   gs = age_h./1000;

    % Nernst
    E = 1.229 - 8.46e-4.*(T - 298.15) + R.*T./(2.*F) .* log(pH2./P0 .* sqrt(pO2./P0));

    % activation (ORR): j0 drops with age
    j0 = exp(P.lnj0 + P.j0f.*gf + P.j0s.*gs) .* exp(-P.Ea.*1e3./R .* (1./T - 1./343.15));
    act = R.*T./(P.alpha.*F) .* (log((j + P.jx)./j0) - P.gamma.*log(pO2./P0));

    % ohmic / membrane: HFR rises linearly with age (imposed from EIS)
    HFR = P.hb(1) + P.hb(2).*age_h./1000 + P.hb(3).*j;
    ohm = j .* HFR .* (opts.RH_c./0.35).^(-P.nu);

    % mass transport: amplitude grows with age
    mt = exp(P.lnc) .* (1 + P.cf.*gf + P.cs.*gs) .* exp(-P.kT.*(T - 343.15)) ...
         .* (j ./ (pO2./P0)).^P.m;

    U = E - act - ohm - mt;
end

% ===================== fitted coefficients =============================
function P = default_params()
    % Fitted by lsqnonlin (sensitivity holds weighted x2); see voltage_model.mat
    P.R = 8.314; P.F = 96485; P.P0 = 101.325; % J/mol/K, C/mol, kPa
    P.Ea = 80; % kJ/mol, ORR activation energy (fixed)
    P.jx = exp(-6); % A/cm^2, crossover-equivalent current (~2.48 mA/cm^2)
    P.tau = 300; % h, fast-degradation time constant
    % voltage-model parameters
    P.lnj0  = -13.6563;   P.j0f = -0.4934;  P.j0s = 0.0801;
    P.alpha =   0.9740;   P.gamma = 0.4081; P.nu  = 0.1985;
    P.lnc   =  -6.9859;   P.m   =  2.0241;
    P.cf    =   0.0897;   P.cs  =  0.2336;  P.kT  = 0.0117;
    % HFR(age,j) = hb(1) + hb(2)*age/1000 + hb(3)*j   [Ohm*cm^2]
    P.hb = [0.0483, 0.0018, 0.0038];
    % thermal model [a, b_flow, c_j2, d_Tin, e_age]
    P.th_T = [-0.3687, 1.6828, 0.2772, -0.0154, -0.1059];
end

function p = psat(t) % saturation pressure of water (Magnus), kPa, t degC
    p = 0.61094 .* exp(17.625 .* t ./ (t + 243.04));
end
function t = dewpoint(pv) % inverse of psat, degC
    l = log(pv ./ 0.61094);
    t = 243.04 .* l ./ (17.625 - l);
end
