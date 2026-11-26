# Framing the optimization problem

With the degradation function established as a plant model, the project turns to a
control scenario: a fuel-cell power generator that charges a battery and holds it near
a target state of charge while the stack ages. The objective is to keep the battery
state of charge (SoC) at 0.92 against a time-varying electrical load, using the fuel
cell's operating setpoints as the manipulated variables.

## Model description

The model has a plant and a controller (Figures 3-6).

**Plant (Figures 3-4).** The fuel cell charges a battery; downstream of the battery the
demand is lumped into a single load signal. The controlled quantity is the battery SoC,
whose dynamics are an energy balance:

    SoC(k+1) = SoC(k) + (Ts / E_cap) * ( P_fc(k) - P_load(k) - P_drain )

where P_fc = N * U * I is the stack power (U and I from the degradation model),
P_load is the measured load, E_cap = V_batt * Q * 3600 is the battery energy capacity,
and Ts the sample time. Because the fuel cell is current-controlled, the controller
sets the current density and the four operating-condition setpoints, not the power
directly. The stack age is advanced by a clock so the degradation model stays active
through the run.

**Why adaptive MPC.** Two features make this a model-predictive control problem rather
than a simple regulator. First, the manipulated variables are bounded (the validated
operating range) and rate-limited (actuator dynamics), and the output must stay in a
band - MPC enforces input and output constraints natively. Second, the plant is
non-stationary: over 2300 h the degradation shifts the current-voltage relationship, so
a fixed internal model would drift. The adaptive MPC re-identifies the plant online - a
recursive ARX estimator updates a low-order model each step from the measured inputs
and SoC, that model is converted to state space and handed to the controller
(Figure 5) - so the predictions follow the aging stack.

**Decision variables, disturbance, output.**

- Manipulated variables (MV):  u = [ current_density, pressure, lambda, temperature, humidity ]   (5)
- Measured disturbance (MD):    d = [ load ]                                                        (1)
- Measured output (OV):         y = [ SoC ]                                                         (1)

The online internal model is

    x(k+1) = A x(k) + B [ u(k); d(k) ],     y(k) = C x(k)

with A, B, C identified online (order matched to the ARX: na = 1, so a single
SoC-like state).

**The receding-horizon optimization.** At each step k the controller solves a
constrained quadratic program over the prediction horizon:

    min over { du(k), ..., du(k+m-1) }

        sum_{i=1..p}  wy^2 ( r - yhat(k+i|k) )^2
      + sum_{i=0..m-1} [ wdu^2 du(k+i)^2  +  wu^2 ( u(k+i) - u_nom )^2 ]

    subject to, along the predicted trajectory:
        x(k+i+1) = A x(k+i) + B [ u(k+i); d(k) ],   yhat = C x
        u_min   <= u(k+i)  <= u_max
        du_min  <= du(k+i) <= du_max
        y_min   <= yhat(k+i) <= y_max

Only the first move u(k) is applied; the horizon recedes one step and the problem is
re-solved. Reference and horizons: r = 0.92, prediction horizon p = 20, control horizon
m = 4, Ts = 1 s.

The output-tracking weight dominates (wy = 10) to hold SoC on target; a small
move-suppression weight (wdu = 0.1) smooths the setpoints; and a small MV weight
(wu = [0, 0.05, 0.05, 0.05, 0.05], toward u_nom) anchors the four weak setpoints near
nominal while leaving current density free.

| Manipulated variable | Min | Max | Rate limit / step |
|---|---|---|---|
| current_density [A/cm^2] | 0.30 | 1.80 | 0.10 |
| pressure [kPaG]          | 50   | 200  | 25   |
| lambda [-]               | 1.6  | 2.1  | 0.25 |
| temperature [degC]       | 55   | 75   | 2.0  |
| humidity [-]             | 0.30 | 0.90 | 0.05 |

Output bound: 0.80 <= SoC <= 1.00.

**Controllability.** SoC is driven almost entirely by current density; pressure,
stoichiometry, temperature and humidity change the cell voltage - and thus the charging
power - only slightly at a given current. In this single-output loop those four are
therefore weakly controllable, which is why they are anchored near their nominal
setpoints by the small MV weight rather than used as active levers. The natural
extension is to add outputs they can regulate - mean cell temperature and a
system-efficiency signal - so the controller trades efficiency against degradation as
the stack ages. That connects back to the net-power framing of the stack-level sweep
(Figure 2): net power, not cell or stack voltage, is the meaningful objective once the
balance-of-plant is charged for.

## Simulate with test load cycle

The closed loop is driven by a repeating step load (Load Consumption, Figure 8) that
switches between a low and a high demand level, so the controller must keep charging the
battery through changing draw.

The response splits into two phases (Figures 7-8). During an initial transient the
adaptive estimator is still learning the plant: the setpoints swing across their full
range and the manipulated variables saturate (current density and lambda ride their
limits), while the online model converges. Once the estimator settles, the loop reaches
steady regulation - SoC is held near the 0.92 target, current density settles to the
low value needed to balance the average load, and pressure, temperature and humidity
converge to a narrow band around their nominal setpoints (about 150 kPaG, 70 degC and
0.35 RH); lambda settles near 2.

Across the run the voltage-loss channel climbs steadily with the age input, confirming
the degradation model remains active and the controller keeps regulating as the stack
ages. The early saturation is the expected cost of starting the online identification
from a seed model; seeding the estimator (initial A, B and a lower initial covariance)
or warm-starting the setpoints shortens that phase.

[Confirm against your scopes: the exact high/low load levels and units, the length of
the transient before it settles, and whether lambda settles at ~2 or nearer the 1.7
command once current density is steady.]
