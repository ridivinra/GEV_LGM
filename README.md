## Code to implement a model, similar to the one introduced in “Approximate Bayesian inference for spatio-temporal flood frequency analysis”

The fitModel.R file will fit the full model from the paper. You might need to adjust and play around with the proposal distribution in the Metropolis-Hastings step.

There is also a more stable model fitted in fitModelNoTrend.R this is the same model except it has no trend in the location parameter of the GEV distribution.

