## Probabilistic jerk finder

To install, you need f2py (part of the standard numpy installation).

1. Clone the repository
2. In the main folder compile the Fortran code with

./make_all

3. The following Jupyter notebooks show how the algorithm works on simple problems: a time series with a single change in slope, and a time series with multiple changes in slope.

```test/test_jerk_finder.ipynb```

```test/test_jerk_finder_multiple.ipynb```





## Function inputs:


 **sigmas**: the perturbations defining the random walk used in the RJMCMC method, in the order sigma_change_value, sigma_move, sigma_birth

**burn_in**: the length of burn in period which is ignored for collecting statistics

**nsample** the length of the Markov chain including the burn_in period

**num_data**: the size of the dataset

**times**: an array of times that define the dataset

**Y**: an array defining the data at the given time (see above)

**delta_Y**: the error in Y, assumed to be the standard deviation of a normal distribution

**Y_min, Y_max**: the bounds of the prior distribution (assumed uniform) for the range of the internal vertices

**times_min, times_max**: the bounds on the model period, which needs to include all the data.

**k_min, k_max**: the bounds on the prior for the number of internal vertices

**discretise_size**: the number of grid points that defines the resolution of some of the outputs (e.g. ensemble mean, median timeseries).

**time_intervals_nbins** the number of bins used to define diagnostics such as the histogram of internal vertices (change points), or the normalised change in slope

**time_intervals_edges**: the edges of the bins used to define diagnostics such as the histogram of internal vertices (change points), or the normalised change in slope.
Note that the size of time_intervals_edges should be time_intervals_nbins + 1

**thin**: the amount of chain thinning. E.g. a value of n means that statistics are only saved for every nth model

**nbins**: the number of bins defining the model density histogram

**credible**: credible interval bounds as percentage e.g. 95. A value of 0 means that credible intervals are not calculated.

**RUNNING_MODE**: 1 indicates that the posterior is returned; otherwise the prior distribution is returned.


## Function outputs

**Acceptance_rates**:  acceptance rates for change value, move (in time), birth, death

**CREDIBLE_SUP**: the upper bound for the credible interval (array size: discretize_size), or an array of zeros if credible interval is not calculated.

**CREDIBLE_INF**: the lower bound for the credible interval (array size: discretize_size), or an array of zeros if credible interval is not calculated.

**ENSEMBLE_AV**: the ensemble average model (array size: discretize_size)

**ENSEMBLE_MEDIAN**: the ensemble median model (array size: discretize_size)

**ENSEMBLE_MODE**: the ensemble modal model (array size: discretize_size)

**CHANGE_POINTS**: a histogram of the timing of internal vertices (array size: time_intervals_nbins). Each model can add at most 1 to each bin; so normalised over all models this is a posterior probability.

**MARGINAL_DENSITY**: a 2D marginal density for the linear ensemble or zeros (if not calculated); (array size: discretise_size,NBINS)

**N_changepoint_hist**: a histogram of the number of internal vertices (1D array, dimensions K_MIN:K_MAX)

**delta_slope** : the sum of the absolute values of any change in slope over all models occuring during a given time bin; normalised by the number of models.
That is, for each time bin: 1/N * \sum_i  abs(change in slope occuring within time bin), where i is over all N models and over a certain temporal bin defined in the input vector time_intervals_edges. (1D array, size: time_intervals_nbins)



