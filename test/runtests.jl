using Test

# Run all test files
@testset "TRPTKiteTurbineSimulator" begin
    include("test_parameters.jl")
    include("test_wind_profile.jl")
    include("test_dynamics.jl")
    include("test_geometry.jl")
end
