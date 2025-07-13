# Spatial-Mark-Resight-Conditional
SMR samplers that update latent individual IDs instead of marginalizing them out. 

We can only marginalize out individual IDs for Poisson count models. By conditioning on a possibly true data set and updating its elements,
or equivalently, updating the individual IDs of each sample, we can do MCMC with any parametric count model. However, when using
the Poisson observation model, the marginal approach in the Spatial-Mark-Resight-Marginal repository is more efficient, so you should 
use that approach. The data format here is the same as SMR Marginal. Currently, I'm adding negative binomial versions of some samplers in the
marginal repository as a model for overdispersion. Luckily, it can still be fit when summing counts over occasions. Zero-truncated Poisson or
negative binomial hurdle models are also good options, but require the full 3D capture history, so are slower to use.

https://github.com/benaug/Spatial-Mark-Resight-Marginal

Latent IDs of each sample are updated one at a time on each MCMC iteration in
the same manner as they are in categorical SMR where each sample can have one or more categorical covariate, like sex or age class. 

https://github.com/benaug/Spatial-Mark-Resight-IDCov

I have another approach working where we update all IDs at each trap at simultaneously, but it is slower than the one at a time approach, at least
with fewer than 1000 samples spread over 81 traps.

So far, this repository contains:
1. Known number of marked individuals (Chandler and Royle 2013, Sollmann et al. 2013), Poisson and negative binomial.

https://www.jstor.org/stable/23566419 https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/12-1256.1

2. Generalized SMR (gSMR) with known number of marked individuals. This includes a marking process to account for different spatial distributions of marked and unmarked individuals (Whittington et al. 2018).generalized SMR (gSMR) with known number of marked individuals. This includes a marking process to account for different spatial distributions of marked and unmarked individuals (Whittington et al. 2018). Poisson and negative binomial observation models.

3. Mb version of 2. Negative binomial observation model only.

https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.12954

