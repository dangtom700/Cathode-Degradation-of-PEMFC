# PEMFC Cathode Degradation Project - Hand-off

MSE 422 course project (Azalea, SFU Mechatronics). Re-analyses the Fensterle et al.
open PEMFC durability dataset, builds a semi-empirical degradation regression model, and
uses it as the plant for an adaptive-MPC control design. Project folder:
`C:\Users\dangs\Desktop\project\fuel_cell_project\Cathode Degradation of PEMFC`.

## Deliverables: TWO separate reports
1. **Analysis report** = `main.mlx` / `main.pdf` (9 pp). The regression / fuel-cell
   analysis. Sections: From the research article, Aim, Methodologies, Analysis
   (degradation over time / spatial / sensitivity / thermal / degradation regression /
   visualizations), Limitations, Further improvement, Conclusion.
   - Closing sections (Limitations / Further improvement / Conclusion) drafted in
     `analysis_report_closing.md` (paragraph style) - paste into main.mlx.
2. **Control-design report** = `control_system_report.md` (standalone). Adaptive MPC that
   charges a battery to a SoC target using the fuel cell. Still needs the real simulation
   figures/results paragraph filled in.
   (`optimization_problem_draft.md` is the older version, now superseded by the control
   report - ignore/delete.)

## The models (numbers already verified this session)

**Degradation voltage model** (`degradation_func.m`; Simulink block version `fcn.m`):
U = E_Nernst - eta_act - eta_ohm - eta_mt. Aging via g_f = 1-exp(-t/300) (fast) and
g_s = t/1000 (slow), acting on j0 (ORR), HFR (ohmic), and MT amplitude.
Fitted: lnj0=-13.66, j0f=-0.493, j0s=0.08, alpha=0.974, gamma=0.408, nu=0.199,
lnc=-6.99, m=2.02, cf=0.090, cs=0.234, kT=0.0117; HFR=0.0483+0.0018*(t/1000)+0.0038*j;
Ea=80 kJ/mol, jx=2.48 mA/cm^2, tau=300 h. Predicts ~34 mV loss at 1.8 A/cm^2 by 2300 h
(ohmic +8, MT +17, kinetics +10). Validity range: p_c 50-200 kPaG, lambda 1.6-2.1,
T 55-75 C, RH 0.30-0.90, j up to ~2.3 A/cm^2.

**Thermal model** (`temperature_model.ipynb` -> `thermal_params.csv`):
T_mean = T_in + a + b*j*6/F_cl + c*j^2 + d*(T_in-70) + e*t/1000; a=-0.369, b=1.683,
c=0.277, d=-0.015, e=-0.106; cross-validated RMSE 0.09 K.

**Stack / battery**: N=10 cells, A=280 cm^2, ~2.6 kW net at BoL, ~6 V stack at full load.
Battery 48 V, 50 Ah ~ 8.64 MJ.

**Key EIS findings** (verified): membrane +8.9% by 1150h, ORR -1.9%->+4.0%, MT ~flat,
total +4.9% at 2300h; spatially outlet-concentrated (+29.8% outlet vs +4.6% inlet total);
pressure sensitivity grows with age; cathode ORR/ECSA is the dominant mechanism.

## Data pipeline (correct run order)
load_data.py -> eis_parser.py -> eis_tier1.py -> clean_load_data.ipynb ->
clean_polarization_data.ipynb -> eis_sensitivity.py -> eis_analysis.ipynb ->
eis_sensitivity.ipynb -> temperature_model.ipynb -> export_matlab.py.
(eis_sensitivity.py MUST follow clean_polarization_data.ipynb - it reads the steady_state
table that notebook builds.) Everything feeds one SQLite DB `cathode_degradation.db`.
MATLAB layer reads `data/holds.csv`, `data/hfr_anchor.csv`, `data/thermal_params.csv`.

## Simulink model (`main.slx`), set up by `controller_setup.m`
Plant: PEMFC Stack (`fcn.m`) -> Dot product (V*I) -> Power Scale -> Battery Status
(outputs dSoC) -> Discrete-Time Integrator -> SoC.
Controller Module: charge controller (SoC->current_density), Air Compressor
(air_flow_demand -> 1/(s+1) lag -> delivered_lambda), Radiator 1/(60s+1), Humidifier
1/(10s+1), each with saturation.
Adaptive MPC: 5 MVs [current_density, pressure, lambda, temperature, humidity], 1 MD
[load], 1 OV [SoC]. Online Recursive Polynomial Model Estimator (ARX, na=1,
nb=nk=ones(1,6)) -> ARX->SS -> one-step Unit Delay on A,B,C,D -> Bus Creator with
A,B,C,D,U,Y,X,DX -> MPC model port. Horizons p=20, m=4. Weights wy=10, wdu=0.1,
wu=[0 0.05 0.05 0.05 0.05]. Ts=1 s (nominal run).

## Gotchas already resolved (do not redo)
- Battery Status outputs the RATE (dSoC), integrated by an external Discrete-Time
  Integrator (IC 0.30, limits [0 1]) - this breaks the algebraic loop. Do NOT update
  persistent state during output.
- `fcn.m` (PEMFC Stack) must NOT have an `arguments` block (Simulink doesn't allow it).
- MPC `model` bus needs all 8 fields; `U` must be a 6x1 2-D matrix (uncheck "Interpret
  vector parameters as 1-D" on the U constant). Name the delayed A/B/C/D signals. Unit
  Delay ICs = seed A0=1, B0=[0.0020 0.0003 0.0002 0.0003 0.0002 -0.0010], C0=1, D0=0.
- ARX order na must equal plant0 state count (=1).
- recursiveARX uses `InitialA`/`InitialB`, NOT `InitialParameters`.
- Estimator Inputs = 6-vector [5 MV + load]; load needs a ZOH at Ts.
- air_flow_demand has NO min-flow floor (floor made lambda saturate at low current).

## Open items / next steps
- **Power Scale** is currently 0.8*3 = 2.4, which caps FC power at ~0.7 kW. Should be
  ~N*eta_DCDC ~ 9 for physical ~2.6 kW. Decide and fix.
- **Age clock** uses 1/86400 (seconds->days, mislabeled "sec>h"). For aging to show over a
  run, use 1/3600 (hours) or time-lapse scale (10 days -> 2300 h via ~1/375).
- **Large-scale run** (`simulate_day_cycle`): Ts=60 s, 10 days, fluctuating load, SoC
  target 1.00. Note target=1.00 sits on the OV upper bound (no headroom) - consider 0.98
  or raise OV.Max; re-identify ARX at Ts=60; keep peak load < FC net capacity.
- **Control report results**: fill the bracketed simulation paragraph with real SoC /
  current / load traces and aging figures (test cycle + 10-day cycle).
- **Weak controllability**: 4 of 5 MVs barely move SoC. Options: drive SoC with
  current_density alone and schedule the others, OR add outputs (cell temperature,
  efficiency) so the other MVs regulate something real.
- Optionally convert control report Limitations/Further work to paragraph style to match
  the analysis report.

## File inventory (project root unless noted)
- degradation_func.m, fcn.m, fc_stack.m, sweep_stack.m, controller_setup.m
- main.mlx, main.pdf, main.slx ; results/ (figures, opt_map.csv, voltage_model.mat)
- data/ : holds.csv, hfr_anchor.csv, thermal_params.csv, thermal_params_boot.csv,
  optimize_operating_point.m
- Notebooks/ : the pipeline .py + .ipynb files above
- cathode_degradation.db (SQLite, large)
- Reports/drafts: control_system_report.md, analysis_report_closing.md,
  optimization_problem_draft.md (superseded)
