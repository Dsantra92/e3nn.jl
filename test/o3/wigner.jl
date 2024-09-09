using e3nn.o3
using Test
using Rotations
using MLUtils: batch

@testset "WignerD" begin
    @testset "basic" begin
        R = rand(RotYXY, 10)
        angles = R .|> Rotations.params
        D = wigner_D.(1, angles)
        @test batch(R - D) .|> abs |> maximum < 1e-8
    end
end
