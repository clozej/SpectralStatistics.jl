
using SpectralStatistics
using Test

@testset "models" begin
    var = [0.545]
    model = Poisson()
    @test level_spacing_pdf(model, var) ≈ [0.5798417833398464]
    @test level_spacing_cdf(model, var) ≈ [0.4201582166601536]
    @test level_spacing_u(model, var) ≈ [0.448952612312622]
    @test number_variance(model, var) ≈ [0.545]
    @test spectral_form_factor(model, var) ≈ [1.0]

    model = GOE()
    @test level_spacing_pdf(model, var) ≈ [0.6779581839643485]
    @test level_spacing_cdf(model, var) ≈ [0.2080704866970233]
    @test level_spacing_u(model, var) ≈ [0.30154207607490247]
    @test number_variance(model, var) ≈ [0.3304590526116955]
    @test spectral_form_factor(model, var) ≈ [0.6882455840426879]

    model = GUE()
    @test level_spacing_pdf(model, var) ≈ [0.6597813324607584]
    @test level_spacing_cdf(model, var) ≈ [0.1401278115132451]
    @test level_spacing_u(model, var) ≈ [0.2442590013361426]
    @test number_variance(model, var) ≈ [0.28956955679634927]
    @test spectral_form_factor(model, var) ≈ [0.545]

    model = GSE()
    @test level_spacing_pdf(model, var) ≈ [0.522349097850527]
    @test level_spacing_cdf(model, var) ≈ [0.06972957445797587]
    @test level_spacing_u(model, var) ≈ [0.17012575578690498]
    @test number_variance(model, var) ≈ [0.2305437559287445]
    @test spectral_form_factor(model, var) ≈ [0.3797911334292492]

end
