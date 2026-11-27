# Optimal Control System Design for a PEMFC Power Generator (MSE 422)

Companion to the cathode-degradation analysis report. That report built a semi-empirical
degradation model of the fuel cell; this report uses it as the plant for a control
design: an adaptive model-predictive controller (MPC) that runs the fuel cell to charge
a battery and hold it near a target state of charge while the stack ages.

## 1. System and scenario

A fuel cell charges a battery that serves a fluctuating electrical load; downstream of
the battery the demand is lumped into a single load signal. The control objective is to
hold the battery state of charge (SoC) at a target (0.92 in the nominal case, 1.00 in
the large-scale case) against the load, using the fuel cell's operating setpoints as the
manipulated variables.

The stack is ~2.6 kW net (10 cells x 280 cm^2). Energy storage is the hydrogen tank plus
a battery buffer; the fuel cell is a charge source, not the store itself. Because a fuel
cell is current-controlled, the controller commands current density and the operating
setpoints, not power directly.

## 2. Plant model

The plant reuses the degradation function from the analysis report (per-cell voltage,
losses, exchange current, HFR and mean cell temperature as functions of current density,
age and operating condition), scaled to the stack and tied to a battery. The controlled
output is SoC, with dynamics given by an energy balance:

    SoC(k+1) = SoC(k) + (Ts / E_cap) * ( P_fc(k) - P_load(k) - P_drain )

where P_fc = N * U * I is the stack power, P_load the measured load, E_cap = V_batt * Q *
3600 the battery energy capacity, and Ts the sample time. The stack age is advanced by a
clock so the degradation model stays active through the run.

## 3. Controller design (adaptive MPC)

Two features make this a model-predictive problem. First, the manipulated variables are
bounded (the validated operating range) and rate-limited (actuator dynamics), and the
output must stay in a band - MPC enforces input and output constraints natively. Second,
the plant is non-stationary: over 2300 h the degradation shifts the current-voltage
relationship, so a fixed internal model would drift. The adaptive MPC re-identifies the
plant online - a recursive ARX estimator updates a low-order model each step from the
measured inputs and SoC, that model is converted to state space and handed to the
controller (one step delayed to break the estimation-control loop), so the predictions
follow the aging stack.

Signals:

- Manipulated variables (MV):  u = [ current_density, pressure, lambda, temperature, humidity ]  (5)
- Measured disturbance (MD):    d = [ load ]                                                       (1)
- Measured output (OV):         y = [ SoC ]                                                        (1)

Online internal model:  x(k+1) = A x(k) + B [ u(k); d(k) ],  y(k) = C x(k), with A, B, C
identified online (ARX order na = 1, a single SoC-like state).

## 4. The optimization problem

At each step k the controller solves a constrained quadratic program over the prediction
horizon:

    min over { du(k), ..., du(k+m-1) }

        sum_{i=1..p}  wy^2 ( r - yhat(k+i|k) )^2
      + sum_{i=0..m-1} [ wdu^2 du(k+i)^2  +  wu^2 ( u(k+i) - u_nom )^2 ]

    subject to, along the predicted trajectory:
        x(k+i+1) = A x(k+i) + B [ u(k+i); d(k) ],   yhat = C x
        u_min   <= u(k+i)  <= u_max
        du_min  <= du(k+i) <= du_max
        y_min   <= yhat(k+i) <= y_max

Only the first move u(k) is applied; the horizon recedes one step and the problem is
re-solved. Reference r = SoC target; prediction horizon p = 20, control horizon m = 4.

Weights: the output-tracking weight dominates (wy = 10) to hold SoC on target; a small
move-suppression weight (wdu = 0.1) smooths the setpoints; a small MV weight
(wu = [0, 0.05, 0.05, 0.05, 0.05], toward u_nom) anchors the weak setpoints near nominal
while leaving current density free.

| Manipulated variable | Min | Max | Rate limit / step |
|---|---|---|---|
| current_density [A/cm^2] | 0.30 | 1.80 | 0.10 |
| pressure [kPaG]          | 50   | 200  | 25   |
| lambda [-]               | 1.6  | 2.1  | 0.25 |
| temperature [degC]       | 55   | 75   | 2.0  |
| humidity [-]             | 0.30 | 0.90 | 0.05 |

Output bound: 0.80 <= SoC <= 1.00.

Controllability. SoC is driven almost entirely by current density; pressure,
stoichiometry, temperature and humidity change the cell voltage only slightly at a given
current, so in this single-output loop they are weakly controllable and are anchored near
nominal rather than used as active levers. The meaningful objective once the
balance-of-plant is charged for is net power (see the stack-level sweep in the analysis
report), which the operating-point optimizer provides as a setpoint schedule.

## 5. Simulation results

Two scenarios were run.

Nominal test cycle (Ts = 1 s). A repeating step load drives the loop. The response has
an initial transient while the online estimator learns the plant - the setpoints swing
and the manipulated variables saturate - then settles to steady regulation: SoC holds
near 0.92, current density settles to the low value that balances the average load, and
pressure, temperature and humidity converge to a narrow band around nominal (about 150
kPaG, 70 degC, 0.35 RH). The voltage-loss channel climbs steadily with the age input,
confirming the degradation model stays active.

Large-scale day cycle (Ts = 60 s, 10 days, target 1.00). A longer run with a more
fluctuating daily load and the battery target raised to full charge, to test
energy-management over life. Note two modelling choices: the age clock must be scaled so
10 real days map to a meaningful fraction of the 2300 h degradation window (otherwise the
stack barely ages over the run); and with the target equal to the SoC upper bound, the
controller chases the constraint with no headroom, so expect current density to ride high
and SoC to hover just below 100%.

[Insert the specific SoC / current-density / load traces and the final aging figures from
the simulation here.]

## 6. Limitations

- Weak controllability: four of the five MVs barely affect SoC, so they do little useful
  work in the single-output loop.
- The identified model is sample-time specific; a model learned at one Ts does not carry
  to another, so the estimator must re-learn when Ts changes.
- Aging is modelled as a function of elapsed time only, not of operating history, so the
  controller cannot yet trade efficiency against degradation rate.
- The voltage model's condition slopes are only trustworthy inside the tested range
  (p_c 50-200 kPaG, lambda 1.6-2.1, T 55-75 degC, RH 0.30-0.90, j up to ~2.3 A/cm^2).
- The startup transient reflects the online identification converging from a seed model.

## 7. Further work

- Add measured outputs the weak MVs can regulate - mean cell temperature and a
  system-efficiency signal - so the full five-input MPC becomes meaningful.
- Identify a degradation-rate model (how operating conditions accelerate aging) to enable
  durability-aware control that trades efficiency against stack life.
- Seed the estimator (initial A, B and lower initial covariance) to shorten the startup
  transient, and validate the closed loop against held-out operating points.

## 8. Conclusion

Using the degradation regression as a plant, an adaptive MPC was designed and simulated
that regulates battery SoC with the fuel cell's setpoints while the stack ages, holding
the target through a fluctuating load. The design demonstrates that the semi-empirical
degradation model is usable as a control-design plant, and sets up the natural next step:
a durability-aware controller that optimises net power and stack life together.
