# legacy/: the original λ/k0 fit

These files reproduce the analysis of the first version of the paper. They are kept for reference only.

- `modifiedPoissonFit.m` and `modelCDF.m`: the weighted-score (WS) fit of λ (units per burst) and k0 (intensity of one unit).
- `ExperimentalAnalysis_WS.m` and `SimulationRunWithAnalysis_WS.m`: the original runners that call the WS fit.

Why the fit is no longer used:
- λ and k0 trade off against each other, and only their product, the mean burst size, is estimated reliably.
- On simulations with known λ, the fit returns λ ≈ 2-3 whatever the true value. The details are in `analysis_outputs/RESULTS.md`, sections 25-31.

Use `GelMetricsReport.m` in the repository root instead. It reports the mean burst size (λ·k0) together with noise-controlled measures and confidence intervals.

The scripts in `checks/` that reproduce the old numbers add this folder to the MATLAB path themselves.
