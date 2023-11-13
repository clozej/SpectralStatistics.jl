
using SpectralStatistics
using Test

@testset "models" begin
    var = [0.545]
    poisson = Poisson()
    @test level_spacing_pdf(poisson, var) ≈ [0.5798417833398464]
    @test level_spacing_cdf(poisson, var) ≈ [0.4201582166601536]
    @test level_spacing_u(poisson, var) ≈ [0.448952612312622]
    @test number_variance(poisson, var) ≈ [0.545]
    @test spectral_form_factor(poisson, var) ≈ [1.0]
end
