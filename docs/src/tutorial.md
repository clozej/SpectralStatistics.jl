# Tutorial
This section will guide you through everything you need to know for getting started with the SpectralStatistics.jl package. To learn more about spectral statistics and quantum chaos in general the books are a good place to start.

## Installation
SpectralStatistics can be installed using the Julia package manager. Simply run:

```@example Install; continued = true  
using Pkg; Pkg.add("SpectralStatistics")
```

This will instal the package to your current environment. To use it run:
```@example Tutorial; continued = true 
using SpectralStatistics
```

in your Julia session.

## Core components
There are three major components of the package:
1. Types that represent spectral data, for example [`UnfoldedSpectrum`](@ref) and [`RealSpectrum`](@ref). These are subtypes of the [`DataSample`](@ref) abstract type.

2. Types that represent analytical models of spectral statistics for example [`Poisson`](@ref) and [`GOE`](@ref). These are subtypes of the [`Model`](@ref) abstract type.

3. Functions that compute the spectral statistics for example [`level_spacing_pdf`](@ref). They can be used either to compute the spectral statistics from a spectrum by imputing the data or evaluate an analytical expression from a model. 

The package also includes convenience codes for fitting the models to data and plotting.

## Preparing and unfolding spectra
Infinite potential wells are usually one of the first example quantum systems one encounters when learning quantum mechanics. We will analyze the spectrum of a two-dimensional well, namely the rectangular quantum billiard. Let us take a rectangle with sides $$a$$ and $$b$$. The energy levels are given by 
```math
E_{n_x, n_y} = \frac{\pi^2}{2}\left(\left(\frac{n_x}{a}\right)^2+\left(\frac{n_y}{b}\right)^2\right),
```
where $$n_x$$ and $$n_y$$ are integer quantum numbers and the particle mass $$m=1$$ and $$\hbar=1$$. We will also fix $$a=1$$ and define the ratio $$\chi=b/a$$. We define a function that returns the first $$N$$ levels:
```@example Tutorial; continued = true
H(nx,ny,a,b) = 0.5*pi^2.0*( (nx/a)^2.0 + (ny/b)^2.0 )

function rectangle_spectrum(N; chi = 1.0*MathConstants.golden)
    M = Int(round(10*chi*sqrt(N)))
    spect = [H(i,j,1.0,chi) for i in 1:M for j in 1:M]
    return sort(spect)[1:N]
end
```
We set the parameter $$\chi$$ to the golden ratio by default. Let us compute the spectrum of the first 10000 levels:
```@example Tutorial; continued = true 
spect = rectangle_spectrum(10000)
```
The SpectralStatistics.jl package uses the `RealSpectrum` container to represent raw spectral data. This is the data we would obtain form an experiment or a numerical simulation. We pass our data to the container:
```@example Tutorial;
real = RealSpectrum(spect)
```
Some spectral statistics may be computed directly from this data. Usually, to compare the universal statistical properties of the spectra of different systems, the spectra must first be normalized in an appropriate way. This process is called spectral unfolding. The density of states varies from system to system. The spectral staircase function number of levels up to some energy
```math
N(E) := \#\{n|E_n<E\}.
```
We are mainly interested in the fluctuations of levels around the local mean that turn out to have universal properties that are related to the dynamics. One way to unfold the spectrum, is to fit a smooth curve $$F(E)$$ to the spectral staircase. This represents the integrated density of states. One of the common ways of unfolding the spectrum is to transform the spectral data by
```math
e_n = F(E_n).
```
The new energies $${e_n}$$ are called the unfolded energies. The function `unfold_spectrum` can be used to unfold a real spectrum. It returns a new spectrum of type `UnfoldedSpectrum`. Polynomials of n-th degrees are commonly used as unfolding functions. We will use a second degree polynomial to unfold the billiard spectrum. We simply call:
 ```@example Tutorial;
unfolded = unfold_spectrum(real, 2)
```
For convenience it is possible to call `unfold_spectrum` with an integer as the second parameter to use a polynomial as the fitting function. In general any julia function (its first parameter must be the energy) may be passed as the second parameter of `unfold_spectrum` and it will be automatically fitted and used as the unfolding function. We have unfolded the spectrum and can now continue with analyzing the spectral statistics.
## Computing spectral statistics
Let us compute one of the most important and commonly used spectral statistics - the nearest-neighbor level spacing distribution. The level spacings are simply the differences between consecutive energies in the spectrum
```math
s_i=e_{i+1}-e_i.
```
If the spectrum is unfolded the mean level spacing will be equal to 1. It is possible to compute the level spacings with the `level_spacing` function. Let us check the unfolding by calling:
 ```@example Tutorial;
using Statistics
mean(level_spacing(unfolded))
```
We see the mean level spacing is indeed approximately equal to 1.
Next we are interested in the probability distribution of finding a certain level spacing in the unfolded spectrum. The `level_spacing_pdf` computes the probability density function (as a histogram) of the level spacing:
 ```@example Tutorial;
s = collect(0.0:0.05:5.0)
p = level_spacing_pdf(unfolded, s)
```
the second parameter must be a `Vector` of evaluation points. We can plot the result by using the excellent [Makie](https://docs.makie.org/stable/) library:
 ```@example Tutorial;
using CairoMakie
f = Figure(resolution = (640,360))
ax = Axis(f[1,1], xlabel=L"s", ylabel=L"P(s)")
lines!(ax,s,p )
save("lsplot.svg", f); nothing # hide
```
![](lsplot.svg)

## Using models
