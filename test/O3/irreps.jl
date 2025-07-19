using E3NN.O3
using Test

@testset "Irreps" begin
    @testset "creation" begin
        Irrep(2, 1)
        ir = Irrep("3o")
        Irrep(ir)
        @test Irrep("10o") == Irrep(10, -1)
        @test Irrep("1y") == Irrep("1o")

        irreps = Irreps(ir)
        Irreps(irreps)
        Irreps([(32, (4, -1))])
        Irreps("11e")

        @test Irreps("16x1e+32x2o") == Irreps([(16, (1, 1)), (32, (2, -1))])
        Irreps(["1e", "2o"])
        Irreps([(16, "3e"), "1e"])
        Irreps([(16, "3e"), "1e", (256, (1, -1))])
    end

    @testset "properties" begin
        irrep = Irrep("3e")
        @test irrep.l == 3
        @test irrep.p == 1
        @test dim(irrep) == 7

        @test Irrep("$(irrep)") == irrep
        ir = Irrep("5o")
        @test ir.l == 5
        @test ir.p == -1

        irreps = Irreps("4x1e + 6x2e + 12x2o")
        @test Irreps(repr(irreps)) == irreps
    end

    @testset "arithmetic" begin
        @test 3 * Irrep("6o") == Irreps("3x6o")
        # We are considering product of irrep as tensor_product
        # products = [ir for ir in Irrep("1o") * Irrep("2e")]
        # @test products == [Irrep("1o"), Irrep("2o"), Irrep("3o")]

        @test Irrep("4o") + Irrep("7e") == Irreps("4o + 7e")
        @test Irreps("1o + 4o") + Irreps("1o + 7e") ==
              Irreps("1o + 4o + 1o + 7e")

        @test 2 * Irreps("2x2e + 4x1o") == Irreps("4x2e + 8x1o")
        @test Irreps("2x2e + 4x1o") * 2 == Irreps("4x2e + 8x1o")
    end

    @testset "empty" begin
        @test Irreps([]) == Irreps("")
        er = Irreps([])
        @test length(er) == 0
        @test dim(er) == 0
        @test ls(er) == []
        @test num_irreps(er) == 0
    end

    @testset "cat" begin
        irreps = Irreps("4x1e + 6x2e + 12x2o") + Irreps("1x1e + 2x2e + 12x4o")
        @test length(irreps) == 6
        @test ls(irreps) == vcat(
            repeat([1], 4),
            repeat([2], 6),
            repeat([2], 12),
            repeat([1], 1),
            repeat([2], 2),
            repeat([4], 12)
        )
        @test lmax(irreps) == 4
        @test num_irreps(irreps) == 4 + 6 + 12 + 1 + 2 + 12
    end

    @testset "contains" begin
        @test Irrep("2e") ∈ Irreps("3x0e + 2x2e + 1x3o")
        @test Irrep("2o") ∉ Irreps("3x0e + 2x2e + 1x3o")
    end

    @testset "errors" begin
        @test_throws ArgumentError Irrep(-1, 1)
        @test_throws ArgumentError Irrep(3, 2)
        @test_throws ArgumentError Irrep("-1e")
        @test_throws ArgumentError Irreps("-1x1e")
        @test_throws ArgumentError Irreps("bla")
    end
end
