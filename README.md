# SpectralStatistics

[![Build Status](https://github.com/clozej/SpectralStatistics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/clozej/SpectralStatistics.jl/actions/workflows/CI.yml?query=branch%3Amain)


Methods for the computation and analysis of the statistics of quantum spectra. 

# Unfolding
The folowing methods for unfolding spectra have been implmented
- fit_integrated_density
- unfold_spectrum
- weyl_law
# List of statistics
Methods for computing the folowing statistics have been implemented:
- mode_fluctuations (mode_fluctuations)
- level spacing distributions (level_spacing, _pdf, _cdf, _u)
- number variance (number_variance)
- spectral form factor (spectral_form_factor)

# List of models
The folowing models have been fully or partialy implemented:
- Poisson
- GOE, GUE, GSE
- Brody
- Berry-Robnik, Berry-Robnik-Brody
- Gamma distributions

# To do list
We would like to include the folowing methods and statistics:
- E(k,l) statistics
- Two point correlation and cluster functions R2, Y2
- spectral rigidity (delta3)
- level spacing ratios and distributions
- complex spectra statistics

We would like to include the folowing models:
- Izrailev level spacing distributions
- non standard RMT ensembles
- More accurate formulas for the RMT ensambles using the Pade aproximant method (Dietz, Haake) 

We would like to implement the integral transformations between the various statistics.

# Dependencies
- Polynomials
- SpecialFunctions
- LsqFit
- CairoMakie for plotting examples