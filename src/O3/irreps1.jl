
"""
    Irreps

Direct sum of irreducible representations of ``O(3)``.

This struct does not contain any data, it is a structure that describes the representation.

# Examples
```jldoctest
# Create a representation of 100 l=0 of even parity and 50 pseudo-vectors.
julia> x = Irreps([(100, (0, 1)), (50, (1, 1))])
100x0e+50x1e

julia> dim(x)
250

# Create a representation of 100 l=0 of even parity and 50 pseudo-vectors.
julia> Irreps("100x0e + 50x1e")
100x0e+50x1e

julia> Irreps("100x0e + 50x1e + 0x2e")
100x0e+50x1e+0x2e

julia> Irreps("100x0e + 50x1e + 0x2e") |> lmax
1

julia> Irrep("2e") in Irreps("0e + 2e")
true

# Empty Irreps
julia> Irreps(), Irreps("")
(, )
```
"""
struct Irreps
    irreps::Tuple{Vararg{MulIrrep}}
end

Irreps(ir::Irrep{T}) where {T <: Integer} = Irreps((MulIrrep(1, ir),))
Irreps(mulir::MulIrrep{T}) where {T <: Integer} = Irreps((mulir,))
function Irreps(mulir::Tuple{Int, Irrep{T}}) where {T <: Integer}
    Irreps((MulIrrep(mulir...),))
end
Irreps() = Irreps{Int}(tuple())

function Base.show(io::IO, irreps::Irreps)
    isempty(irreps.irreps) && return print(io, "")
    print(io, join(string.(irreps.irreps), " + "))
end

Base.iterate(irreps::Irreps, args...) = iterate(irreps.irreps, args...)
Base.getindex(irreps::Irreps, i) = getindex(irreps.irreps, i)
Base.length(irreps::Irreps) = length(irreps.irreps)
Base.in(ir::Irrep, irreps::Irreps) = any(mulir -> mulir.ir == ir, irreps.irreps)

function Base.:+(a::Irreps, b::Irreps)
    return regroup(Irreps((a.irreps..., b.irreps...)))
end

Base.:+(a::Irrep, b::Irrep) = Irreps(a) + Irreps(b)

function Base.:*(n::Integer, irreps::Irreps)
    @assert n >= 0 "Multiplicity must be non-negative."
    return Irreps(tuple([(n * mul, ir) for (mul, ir) in irreps]...))
end

function Base.:*(n::Integer, ir::Irrep)
    @assert n >= 0 "Multiplicity must be non-negative."
    return Irreps(((n, ir),))
end

Base.:*(ir::Irrep, n::Integer) = n * ir

function Base.div(irreps::Irreps, n::Integer)
    Irreps(tuple([(div(mul, n), ir) for (mul, ir) in irreps]...))
end
function Base.repeat(irreps::Irreps, n::Integer)
    Irreps(tuple(repeat(collect(irreps.irreps), n)...))
end

function Base.show(io::IO, irreps::Irreps)
    print(io, join(["$(mul)x$(ir)" for (mul, ir) in irreps.irreps], " + "))
end

dim(irreps::Irreps) = sum(mul * dim(ir) for (mul, ir) in irreps; init = 0)
num_irreps(irreps::Irreps) = sum(first, irreps.irreps; init = 0)
ls(irreps::Irreps) = [ir.l for (mul, ir) in irreps for _ in 1:mul]
function lmax(irreps::Irreps)
    isempty(irreps.irreps) ? -1 : maximum(ir.l for (_, ir) in irreps)
end
function count(ir::Irrep, irreps::Irreps)
    sum(mul for (mul, irr) in irreps if irr == ir; init = 0)
end
function mul_gcd(irreps::Irreps)
    isempty(irreps.irreps) ? 0 : gcd(first(x) for x in irreps)
end
isscalar(irreps::Irreps) = all(ir.l == 0 && ir.p == 1 for (_, ir) in irreps)

function remove_zero_multiplicities(irreps::Irreps)
    Irreps(tuple([(mul, ir) for (mul, ir) in irreps if mul > 0]...))
end

function unify(irreps::Irreps)
    out = Vector{Tuple{Int, Irrep}}()
    for (mul, ir) in irreps
        if !isempty(out) && last(out)[2] == ir
            out[end] = (last(out)[1] + mul, ir)
        else
            push!(out, (mul, ir))
        end
    end
    return Irreps(tuple(out...))
end

simplify(irreps::Irreps) = unify(remove_zero_multiplicities(irreps))

function regroup(irreps::Irreps)
    counts = Dict{Irrep, Int}()
    for (mul, ir) in irreps
        counts[ir] = get(counts, ir, 0) + mul
    end
    sorted_irreps = sort(collect(keys(counts)))
    return Irreps(tuple([(counts[ir], ir)
                         for ir in sorted_irreps if counts[ir] > 0]...))
end

function _sort_with_p_info(irreps::Irreps)
    v = collect(irreps.irreps)
    p = sortperm(v, by = x -> x[2])
    inv_p = invperm(p)
    sorted_irreps = Irreps(tuple(v[p]...))
    return (irreps = sorted_irreps, p = p, inv = inv_p)
end

function Base.sort(irreps::Irreps)
    return _sort_with_p_info(irreps)[irreps]
end
