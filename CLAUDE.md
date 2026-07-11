# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo does

MATLAB code for estimating Kilauea's pressure/volume history during the 2018 collapse from a joint inversion of GPS, tilt (SDH), and InSAR (ascending + descending) data. The forward model represents two pressurized spheroids — the Halemaʻumaʻu (HMM) and South Caldera (SC) chambers — via Yang's solution (`spheroid.m`). MCMC (Metropolis-Hastings, `mcmc.m`) and pattern-search L-curve sweeps (`generate_l_curve_point.m`) jointly invert chamber geometry and pressure/volume change parameters, after which `TimeDependentLSQtilt.m` solves the time-dependent pressure history per timestep.

## Common workflows (in MATLAB)

These scripts are run interactively in MATLAB from the repo root — there is no build/test/lint system.

- **Build priors (run once, before any MCMC):** `prior_prob` — reads PNGs in `Data/post_im/`, kernel-density-fits Normal/Lognormal/Gamma/Uniform distributions, and writes `Data/paramDists.mat`. All downstream MCMC/L-curve code loads this file.
- **Build InSAR covariance matrices:** scripts in `Make_cov_insar/` (`Prep_InSAR_cov.m`) — produces `Data/asc_cov.mat` and `Data/desc_cov.mat`.
- **Main driver — geometry inversion + pressure history:** `optimize_geo_weighted_MCMC` — top-level script that loads GPS/tilt/InSAR, runs MCMC or an L-curve sweep, then calls `TimeDependentLSQtilt` and the plotting functions. Toggle behavior via the flags near line ~303: `solveweights`, `runMCMC`, `run_L_curve`, `l_curve_type` (`"prior"` or `"gps"`).
- **Plotting:** scripts in `Plotting functions/` (`plotHists`, `plot_insar`, `plot_insar_new`, `makeplots`, `plotCorr`, `plotLcurve`, `l_curve_test`) — invoked from the main driver after a run completes.

`addpath('./Plotting functions/')` is added inside the main driver; the `Make_cov_insar/` folder is run separately.

## Architecture — pieces that span multiple files

### Parameter vector conventions (load-bearing — read before editing)

The MCMC parameter vector `opt_params` is 16 elements, in this fixed order (see `optimize_geo_weighted_MCMC.m` ~line 275 and `get_full_m.m`):

```
[dvHMM_insar, dvHMM_gps, volHMM, xHMM, yHMM, dHMM, alphaHMM,
 xSC,        ySC,        dSC,    alphaSC, dipSC, strikeSC,
 dvSC_insar, dvSC_gps,   volSC]
```

`get_full_m(prior, opt_params, forward, pressure_type)` is the canonical converter between this MCMC layout and the 1×16 `m` vector that `spheroid.m` expects (`[vert_sd, horiz_sd, dip, strike, x, y, z, strength] × 2`). Two distinct chamber pressure/volume changes (`_gps` vs `_insar`) are carried through the chain because GPS and InSAR sample different time windows of the eruption — `pressure_type` selects which one to use when building the forward model. **Any change to the parameter ordering must be mirrored in `get_full_m.m`, `mcmc.m` (`priorNames` lookup), `create_MCMC_data.m`, `prior_prob.m` (`paramDists` field names), `check_intersection.m`, and the `paramNames`/`lb`/`ub` arrays in `optimize_geo_weighted_MCMC.m`.**

`taiyi_parameters` (16-element vector from Wang et al. 2021) is the standard warm-start / prior geometry — it appears literally in both `mcmc.m` and `optimize_geo_weighted_MCMC.m`. If it is edited, update both.

### Data pipeline

1. `loadinsar.m`, `coord.m`, `llh2local.m`, `xyz2llh.m` — geographic + coordinate utilities; the local origin used throughout is `[-155.2784, 19.4073]` (Kilauea summit, in km then scaled to m).
2. `create_MCMC_data.m` — forward model: builds the stacked data vector `[GPS_data; insar_data]` from an `opt_params` proposal. The MCMC and L-curve cost functions both call this; the `length(varargin{3})*3` GPS-block-length convention in `mcmc.m` assumes this layout.
3. `mcmc.m` — Metropolis-Hastings sampler. Two modes via `solveweights`: fixed weights (uses `gps_weight`, `insar_weight`, `prior_weight` directly) or hierarchical weights (last two params of `x` are `log_gamma_gps`, `log_gamma_insar` and the chain samples them). Proposals are rejected if either (a) they exit `xbnds` or (b) `check_intersection` says the two source bounding spheres overlap — adding/removing this constraint changes the geometry that the chain can explore.
4. `generate_l_curve_point.m` — wraps `patternsearch` (Global Optimization Toolbox) for L-curve sweeps. Faster than running an MCMC per weight value and is the current default (see commit `b59d5cd`).
5. `TimeDependentLSQtilt.m` — once geometry is fixed, solves a single large weighted LSQ for `[dpHMM(t); dpSC(t); station_offsets]` across all timesteps, with a random-walk regularizer (`rwsigma` from `GetRandomWalk`) tying neighboring `dp` values.
6. `GetErrors.m` — Monte Carlo error propagation: draws `N_draws` geometries from the posterior, adds `N_noise` realizations of synthetic GPS noise, re-runs the LSQ to bracket `dp_low/dp_high` and `u_low/u_high`.

### Forward model

`spheroid.m` (Yang et al.) — internal/surface deformation from a dipping pressurized spheroid in a half-space. `strength_type` argument toggles between pressure-change and volume-change input; the MCMC path uses `'volume'`, while the geometry-only LSQ path used `'pressure'`. `spheroid_pFromV.m` converts between the two given Poisson ratio and shear modulus (`0.25`, `3.08e9` are hard-coded throughout). `creategreens.m` / `createtiltgreens.m` evaluate the displacement and tilt Green's functions at the GPS and tilt stations for a frozen geometry.

### Data layout

- `Data/` — input/output `.mat` files, posterior images, processed GPS/InSAR/tilt time series. Gitignored.
- `Make_cov_insar/` — InSAR covariance estimation (variogram fits, ASCII InSAR text inputs).
- `Figures/`, `PaperFigs/` — output PDFs/PNGs from plotting scripts. Gitignored.
- `Elizabeth insar scripts/` — external InSAR preprocessing scripts, not part of the main pipeline. Gitignored.
- `MCMC Test/` — sandbox / older MCMC experiments.

`Data/`, `Figures/`, `PaperFigs/`, `Elizabeth insar scripts/`, and all `*.mat` files are gitignored — large generated artifacts live alongside the code but are not checked in.

### Dependencies

- MATLAB with Statistics and Machine Learning Toolbox (`pdf`, `ksdensity`, `mvnrnd`, etc.).
- Global Optimization Toolbox for `patternsearch` (used by `generate_l_curve_point`).
- Parallel Computing Toolbox is optional — `UseParallel` is set true in `patternsearch` options.

## Workflow gotchas

- `prior_prob.m` must be run (and `Data/paramDists.mat` regenerated) before any MCMC, or `mcmc.m` will fail on `load Data/paramDists.mat`.
- Several `priorNames` entries in `prior_prob.m` reuse the same source image (e.g. `alpha_sc.png` is listed for `dpHMM_insar`, `dvHMM_gps`, `volSC`, etc.) — these are deliberately broad fallbacks, not real priors for those parameters. Don't "fix" them without checking what each name is actually used for downstream.
- MCMC outputs are saved as `.mat` files in `Data/` with descriptive names like `MCMC_final_2e6_dv_gps26_insar016_nointersect.mat`; the main driver loads one by name when `runMCMC = false`. Update the filename when iterating.
- The `nanstatbeginning`/`nanstatend` masks drop GPS stations CRIM, UWEV, BYRL (and CALS) from the inversion in multiple places — keep them consistent if changing the station list.
