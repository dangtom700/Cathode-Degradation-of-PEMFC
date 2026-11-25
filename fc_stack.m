function s = fc_stack(j, age_h, x, P)
    % Stack voltage/power + compressor & pump at operating point x, nominal setpoint.
    pc = x(1); lam = x(2); Tin = x(3); RH = x(4);
    
    % per-cell voltage from the degradation model (computes T_cell internally)
    out = degradation_func(j, age_h, 'p_c', pc, 'lambda_c', lam, 'T_in', Tin, 'RH_c', RH, ...
        'lambda_a', P.lam_a, 'RH_a', P.RH_a, 'dpa', P.dpa, 'F_cl', P.F_cl);
    U = out.U;  T_cell = out.T_cell;
    
    I = j * P.A;                                              % A
    nO2 = P.N * I / (4*P.F);  nair = lam * nO2 / 0.2095;      % mol/s (dry air)
    mair = nair * 0.02897;                                    % kg/s
    Fn = nair * 22.414 * 60;                                  % NL/min
    
    Pout = pc + P.P0;  dp = 19;
    for it = 1:3                                              % stack drop at mean channel pressure
        Qv = Fn * P.P0 / (Pout + dp/2) * (Tin + 5 + 273.15) / 273.15;
        dp = P.dp_coef(1)*Qv + P.dp_coef(2)*Qv^2;
    end
    P_cmp_out = Pout + dp + P.dp_hum;
    
    s.P_comp  = mair * 1005 * P.T_amb / (P.eta_is * P.eta_mot) * ((P_cmp_out/P.P0)^(0.4/1.4) - 1);
    s.P_pump  = P.dp_cl * (P.F_cl/6)^2 * 1e3 * P.F_cl/60000 / P.eta_pump;
    s.P_stack = P.N * I * U;
    s.P_net   = s.P_stack - s.P_comp - s.P_pump;
    s.V_stack = P.N * U;
    s.U = U;  s.T_cell = T_cell;
end