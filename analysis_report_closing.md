# Closing sections for main.pdf (analysis report)

## Limitations

- Single-measurement EIS slices: most sensitivity spectra are one measurement per
  condition, so the ratios are indicative and carry no error bars.
- Two-stack BoL comparison: the BoL variation spectra come from a fresh 8-cell stack, not
  the 10-cell durability stack, so BoL-vs-aged comparisons use within-stack ratios rather
  than absolute values.
- The HFR age slope is imposed from the Tier-1 EIS, so the voltage model presents the
  ohmic degradation rather than independently fitting it.
- Aging is captured as a function of elapsed time only, not of operating history, so the
  model describes how the stack has aged, not how operating conditions change the aging
  rate.
- The condition slopes (gamma, nu, kT, m) are only trustworthy inside the tested range:
  p_c 50-200 kPaG, lambda 1.6-2.1, T 55-75 degC, RH 0.30-0.90, j up to ~2.3 A/cm^2.
- The small negative thermal aging term (-0.1 K per 1000 h) is most likely sensor drift,
  not physical, and is kept only to keep the residuals clean.

## Further improvement

- Add repeats and error bars to the sensitivity EIS so the operating-condition slopes
  carry uncertainty.
- Identify a degradation-rate model - how pressure, temperature, stoichiometry and
  humidity accelerate aging - so the regression can support durability-aware operation.
- Extend the Tier-2 per-cell circuit fits and couple the current-density map features
  into the voltage model.
- Validate the voltage model against held-out ages and operating points.
- Apply the model: the optimal control system design built on this plant is covered in a
  separate report.

## Conclusion

This project re-analysed the Fensterle et al. durability dataset and condensed it into a
usable model. The cleaning pipeline rebuilt consistent tables from the raw logs,
polarization curves and EIS spectra; the EIS analysis reproduced the paper's findings on
when and where the cathode degrades (ORR charge-transfer dominated, concentrated at the
outlet, sharpest at low pressure, with the membrane and mass transport nearly flat); and
a semi-empirical voltage model folded all of this into explicit aging terms, predicting
about 34 mV of loss at 1.8 A/cm^2 by 2300 h. Treating age as a continuous parameter
rather than a BoL/aged label gives a physics-informed plant suitable for operating-point
optimization and control, which the companion control-design report builds on.
