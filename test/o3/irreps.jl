using e3nn.o3
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

        # This is a weird way to do this,
        # Iterators.take is a better way to do this
        # Just wanted to test if the iterator works infinitely
        iter = Irrep
        for x in range(0, 500)
            irrep, iter = Iterators.peel(iter)
            @test irrep.l == x // 2 |> trunc
            @test irrep.p in (-1, 1)
            @test dim(irrep) == (2 * trunc(x // 2) + 1)
        end

        irreps = Irreps("4x1e + 6x2e + 12x2o")
        @test Irreps(repr(irreps)) == irreps
    end

    @testset "arithmetic" begin
        @test 3 * Irrep("6o") == Irreps("3x6o")
        products = [ir for ir in Irrep("1o") * Irrep("2e")]
        @test products == [Irrep("1o"), Irrep("2o"), Irrep("3o")]

        @test Irrep("4o") + Irrep("7e") == Irreps("4o + 7e")
        @test Irreps("1o + 4o") + Irreps("1o + 7e") == Irreps("1o + 4o + 1o + 7e")

        @test 2 * Irreps("2x2e + 4x1o") == Irreps("2x2e + 4x1o + 2x2e + 4x1o")
        @test Irreps("2x2e + 4x1o") * 2 == Irreps("2x2e + 4x1o + 2x2e + 4x1o")
        @test Irreps("2x2e + 4x1o") * 2 |> regroup == Irreps("8x1o + 4x2e") # note the ordering
    end

    @testset "empty" begin
        @test Irreps([]) == Irreps("")
        er = Irreps([])
        @test length(er) == 0
        @test dim(er) == 0
        @test ls(er) == []
        @test num_irreps(er) == 0
    end

    @testset "getitem" begin
        irreps = Irreps("16x1e + 3e + 2e + 5o")
        @test irreps[1] == MulIrrep(16, Irrep("1e"))
        @test irreps[4] == MulIrrep(1, Irrep("5o"))
        @test irreps[length(irreps)] == MulIrrep(1, Irrep("5o"))
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

    @testset "ordering" begin
        n_test = 100

        last = nothing
        for (i, irrep) in enumerate(Iterators.take(Irrep, n_test))
            if !isnothing(last)
                @test last < irrep
            end
            if i == n_test
                break
            end
            last = irrep
        end
    end

    @testset "empty" begin
        @test Irreps() == Irreps("")
        @test Irreps("") == Irreps([])
        @test length(o3.Irreps()) == 0
        @test dim(Irreps()) == 0
        @test ls(Irreps()) == []
        @test num_irreps(Irreps()) == 0
    end

    @testset "getitem" begin
        irreps = Irreps("16x1e + 3e + 2e + 5o")
        @test irreps[1] == (16, Irrep("1e"))
        @test irreps[4] == (1, Irrep("5o"))
        @test irreps[end] == (1, Irrep("5o"))

        sliced = irreps[3:end]
        @test sliced isa Irreps
        @test sliced == Irreps("2e + 5o")
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
