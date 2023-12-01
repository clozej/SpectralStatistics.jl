# Tutorial
This section will guide you through everything you need to know for getting started with the SpectralStatistics.jl package.

!!! info "Quantum chaos"
    To learn more about spectral statistics and quantum chaos in general the books *Quantum Chaos: An Introduction* by H.-J. Stöckmann  and *Quantum Signatures of Chaos* by F. Haake are a good place to start.

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
We set the parameter $$\chi$$ to the golden ratio by default. Let us compute the spectrum of the first 50000 levels:
```@example Tutorial; continued = true 
spect = rectangle_spectrum(50000)
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
lines!(ax,s,p)
save("lsplot.svg", f); nothing # hide
```
![](lsplot.svg)

## Using models
We will now compare the results with some analytical models. The library features many of the most used ones. For a complete list see the [API](@ref) section.
Let us compare the level spacings with the Poisson model. To initialize a model we call its constructor:
 ```@example Tutorial;
poisson = Poisson()
```
We can check the characteristics of the model in the documentation [`Poisson`](@ref), where we read that it describes the spectral statistics of an integrable system.  
This is a parameter-less model so the constructor takes no arguments. The level spacing probability density function for the Poisson model is an exponential decay
```math
P(s)=\exp(-s).
```
We can evaluate the analytical expression by calling the `level_spacing_pdf` on the model in the same way as we did for the spectral data:
 ```@example Tutorial;
p_poisson = level_spacing_pdf(poisson, s)
```
Let us add the result to the plot:
 ```@example Tutorial;
lines!(ax,s,p_poisson)
save("lsplot2.svg", f); nothing # hide
```
![](lsplot2.svg)

We see the model curve fits the data nicely and we can conclude the rectangular billiard belongs to the class of integrable systems. If this was a real application we might be satisfied with the analysis. However, just to demonstrate how to fit more complicated models we will try describing the data with another. Many models are actually parametric families, where the parameters are given as fields of the Julia type that codes the model. Let's use the [`Brody`](@ref) model. From the documentation we see this model interpolates between the `Poisson` and `GOE` models as the parameter `:beta` goes form 0.0 to 1.0. We initialize the model by calling:

 ```@example Tutorial;
brody = Brody()
```
Here we initialized the model with the default value of the parameter. We can see from the output this is equal to 1.0. To initialize the model with a different value of the parameter we simply call the constructor with the selected value as an argument: 

 ```@example Tutorial;
brody = Brody(0.2)
```

Some models might posses several parameters. One example is the [`BerryRobnikBrody`](@ref) model. To check all the names of the parameters we can use the `fieldnames` function on the associated type:
 ```@example Tutorial;
fieldnames(BerryRobnikBrody)
```
We see the parameters of this model are called `rho` and `beta`. To make sure we do not mix up the order of the parameters we may also initialize the model with keyword arguments corresponding to the names of the parameters:
 ```@example Tutorial;
brb = BerryRobnikBrody(beta=0.1,rho=0.5)
```
Note from the output that the parameters are initialized correctly regardless of the input order.
Usually we want to fit the model to some data in order to find the optimal values of the parameters. We can do this with the `fit_model` function. This will return a new instance of the model with the adjusted fitting parameter. Since, we have the level spacings already stored we can call:
 ```@example Tutorial;
brody = fit_model(brody, level_spacing_pdf, s, p)
```
Here we chose to overwrite the variable `brody` with the a new instance with an adjusted parameter value. The second parameter is the statistic we wish to fit.
We may also call:
 ```@example Tutorial;
brody = fit_model(unfolded, brody, level_spacing_pdf)
```
to compute the statistic directly from the spectrum if we did not do so beforehand. We see that the fitting parameter is almost 0 and it agrees with the Poisson model.

To learn more about all the features of the library continue by exploring the [API](@ref) section.