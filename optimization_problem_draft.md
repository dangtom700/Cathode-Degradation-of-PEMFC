# Framing the optimization problem

With the degradation function established as a plant model, the project turns to a
control scenario: a fuel-cell power generator that charges a battery and holds it near
a target state of charge while the stack ages. The objective is to keep the battery
state of charge (SoC) at 0.92 against a time-varying electrical load, using the fuel
cell's operating setpoints as the manipulated variables. The full model is shown in
Figures 3-6 and the closed-loop response in Figures 7-8.

## Control scenario and plant

The fuel cell charges a battery; downstream of the battery the demand is lumped into a
single load signal. The controlled quantity is the battery SoC, whose dynamics are an
energy balance:

    SoC(k+1) = SoC(k) + (Ts / E_cap) * ( P_fc(k) - P_load(k) - P_drain )

where P_fc = N * U * I is the stack power (U and I from the degradation model),
P_load is the measured load, E_cap = V_batt * Q * 3600 is the battery energy capacity,
and Ts the sample time. Because the fuel cell is current-controlled, the controller
sets the current density and the four operating-condition setpoints, not the power
directly.

## Why adaptive MPC

Two features make this a model-predictive control problem rather than a simple
regulator. First, the manipulated variables are bounded (the validated operating
range) and rate-limited (actuator dynamics), and the output must stay in a band - MPC
enforces input and output constraints natively. Second, the plant is non-stationary:
over 2300 h the degradation shifts the current-voltage relationship, so a fixed
internal model would drift. The adaptive MPC re-identifies the plant online - a
recursive ARX estimator updates a low-order model each step from the measured inputs
and SoC, that model is converted to state space and handed to the controller
(Figure 5) - so the predictions follow the aging stack.

## Decision variables, disturbance, output

- Manipulated variables (MV):  u = [ current_density, pressure, lambda, temperature, humidity ]   (5)
- Measured disturbance (MD):    d = [ load ]                                                        (1)
- Measured output (OV):         y = [ SoC ]                                                         (1)

The online internal model is

    x(k+1) = A x(k) + B [ u(k); d(k) ],     y(k) = C x(k)

with A, B, C identified online (order matched to the ARX: na = 1, so a single
SoC-like state).

## The receding-horizon optimization

At each step k the controller solves a constrained quadratic program over the
prediction horizon:

    min over { du(k), ..., du(k+m-1) }

        sum_{i=1..p}  wy^2 ( r - yhat(k+i|k) )^2
      + sum_{i=0..m-1} [ wdu^2 du(k+i)^2  +  wu^2 ( u(k+i) - u_nom )^2 ]

    subject to, along the predicted trajectory:
        x(k+i+1) = A x(k+i) + B [ u(k+i); d(k) ],   yhat = C x
        u_min   <= u(k+i)  <= u_max
        du_min  <= du(k+i) <= du_max
        y_min   <= yhat(k+i) <= y_max

Only the first move u(k) is applied; the horizon recedes one step and the problem is
re-solved (receding horizon).

Reference and horizons:  r = 0.92,  prediction horizon p = 20,  control horizon m = 4,
Ts = 1 s.

Weights. The output-tracking weight dominates (wy = 10) to hold SoC on target; a small
move-suppression weight (wdu = 0.1) smooths the setpoints; and a small MV weight
(wu = [0, 0.05, 0.05, 0.05, 0.05], toward u_nom) anchors the four weak setpoints near
nominal while leaving current density free.

Constraints (physical units):

| Manipulated variable | Min | Max | Rate limit / step |
|---|---|---|---|
| current_density [A/cm^2] | 0.30 | 1.80 | 0.10 |
| pressure [kPaG]          | 50   | 200  | 25   |
| lambda [-]               | 1.6  | 2.1  | 0.25 |
| temperature [degC]       | 55   | 75   | 2.0  |
| humidity [-]             | 0.30 | 0.90 | 0.05 |

Output bound: 0.80 <= SoC <= 1.00.

## Controllability note

SoC is driven almost entirely by current density; pressure, stoichiometry, temperature
and humidity change the cell voltage - and thus the charging power - only slightly at a
given current. In this single-output loop those four are therefore weakly controllable,
which is why they are anchored near their nominal setpoints by the small MV weight
rather than used as active levers. The natural extension is to add outputs they can
regulate - mean cell temperature and a system-efficiency signal - so the controller
trades efficiency against degradation as the stack ages. That connects back to the
net-power framing: the stack-level sweep (Figure 2) showed net power, not cell or stack
voltage, is the meaningful objective once the balance-of-plant is charged for, and the
operating-point optimizer (optimize_operating_point.m) produces the setpoint schedule
those four MVs would track.

## Closed-loop behaviour

Figures 7-8 show the controller holding SoC near 0.92 against a square-wave load
(Load Consumption, Figure 8): current density rises when the load draws the battery
down and falls as SoC recovers, while pressure, temperature and humidity stay within a
narrow band around their nominal values. The voltage-loss channel climbs steadily with
the age input, confirming the degradation model remains active across the run.
