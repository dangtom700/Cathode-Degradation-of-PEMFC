# Closing sections for main.pdf (analysis report) - paragraph style

## Limitations

The regression is a semi-empirical model fit to a single open dataset, and that limits
how strongly it can be validated. The dataset covers one 2300-hour run under cyclic
loading with essentially no repeated measurements, so there is no independent data to
confirm that the fitted regression generalises rather than simply interpolating this one
history. The original study used the data to establish which degradation mechanisms are
significant and where they act, but placed little emphasis on when they develop over
time - yet the time axis is exactly what this project's regression depends on, so the
aging terms are supported by the trend of a single run rather than by repeated,
time-resolved validation. Finally, the model leans heavily on the Nernst equation and
other empirical sub-models to represent the cell's behaviour over long operation, so its
accuracy is tied to assumptions - a fixed activation energy, an imposed HFR trend, and
validity only within the tested operating range - that were not independently verified
here.

## Further improvement

The clearest improvement is validation data. Repeated runs, or measurements from a
second stack aged under the same load cycle, would show whether the time-dependent terms
hold up or merely fit this one history; intermediate-age checkpoints matter most, since
the aging factors are what the dataset constrains least. Beyond validation, the model
only describes how the stack has aged as a function of elapsed time, so a dedicated
degradation-rate model - how pressure, temperature, stoichiometry and humidity accelerate
or slow the aging - would let the regression predict rather than replay, and would
support operating the stack for longer life. Reducing the reliance on the empirical
sub-models, for example by fitting the ohmic term directly instead of imposing the HFR
trend, and testing the fit against held-out operating points, would strengthen it
further. The model is then ready to be applied: the optimal control system design built
on this plant is covered in a separate report.

## Conclusion

This project re-analysed the Fensterle et al. durability dataset and condensed it into a
usable model. The cleaning pipeline rebuilt consistent tables from the raw load-cycle
log, polarization curves and EIS spectra; the EIS analysis reproduced the study's
findings on when and where the cathode degrades - charge-transfer (ORR) dominated,
concentrated at the cathode outlet, sharpest at low pressure, with the membrane and
mass transport nearly unchanged - and a semi-empirical voltage model folded these
effects into explicit aging terms, predicting roughly 34 mV of loss at 1.8 A/cm^2 by
2300 hours. By treating age as a continuous parameter rather than a beginning-of-life or
aged label, the model becomes a physics-informed plant that can be swept, optimised and
controlled, even if, as noted, its time-dependent behaviour rests on a single unrepeated
run. That plant is the foundation for the optimal control system design developed in the
companion report.
