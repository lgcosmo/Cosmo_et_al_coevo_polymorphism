using DrWatson, Test
@quickactivate "Cosmo_et_al_coevo_polymorphism"

# Here you include files using `srcdir`
# include(srcdir("file.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "Cosmo_et_al_coevo_polymorphism tests" begin
    @test 1 == 1
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
