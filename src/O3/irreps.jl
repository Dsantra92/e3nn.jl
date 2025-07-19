import Base.Iterators

"""
Irrreducible representations of ``O(3)``.

This struct does not contain any data; it is a structure that describes the representation.
It is typically used as an argument in other parts of the library to define the input and output representations of
functions.

# Fields:
- `l::Int`: non-negative integer, the degree of the representation, ``l = 0, 1, \\dots``
- `p::Int`: the parity of the representation, either 1 (even) or -1 (odd)

The degree ``l`` corresponds to the angular momentum quantum number and defines the dimension of the irrep, with ``l = 0`` representing a scalar, ``l = 1`` a vector, and higher ``l`` values corresponding to tensor-like objects.
"""
struct Irrep{T <: Integer}
    l::T
    p::T

    @doc"""
         Irrep(l::Integer, p::Integer)

     Instantiate a new [`O3.Irrep`](@ref) object.

     # Examples:
     Create a scalar representation (``l=0``) of even parity:
     ```jldoctest
     julia> Irrep(0, 1)
     0e
     ```

     Create a pseudotensor representation (``l=2``) of odd parity:
     ```jldoctest
     julia> Irrep(2, -1)
     2o
     ```
     """
    function Irrep(l::Integer, p::Integer)
        (l >= 0) ||
            throw(ArgumentError("l must be zero or positive integer, got $l"))
        (p in (-1, 1)) ||
            throw(ArgumentError("parity(p) must be one of (-1, 1) got $p"))
        promoted_l, promoted_p = promote(l, p)
        return new{typeof(promoted_l)}(promoted_l, promoted_p)
    end
end

"""
    Irrep(ir::T) where {T <: AbstractString}

Instantiate a new [`O3.Irrep`](@ref) object from it's string representation.

The string representation should be of the form `l` followed by "e" or "o" for even or odd parity, respectively.
There can also be a "y" at the end, which is used to represent the parity of the spherical harmonics.

- "e" for even parity, translates to `p`=1
- "o" for odd parity, translates to `p`=-1
- "y" for the parity of the spherical harmonics, translates to `p`=``(-1)^l``

# Examples
Create a vector representation (``l=1``) of the parity of the spherical harmonics (``-1^l``) gives odd parity):
```jldoctest
julia> Irrep("1y")
1o
```
"""
function Irrep(ir::AbstractString)::Irrep{Int}
    name = strip(ir)
    try
        l = parse(Int, name[1:(end - 1)])
        (l >= 0) ||
            throw(ArgumentError("l must be zero or positive integer, got $l"))
        p = Dict('e' => 1, 'o' => -1, 'y' => (-1)^l)[name[end]]
        return Irrep(l, p)
    catch
        throw(ArgumentError("Cannot convert string $name to an Irrep"))
    end
end

"""
    Irrep(ir::Tuple{Integer, Integer})

Instantiate a new [`O3.Irrep`](@ref) object from a tuple of (`l`, `p`).
"""
Irrep(ir::Tuple{Integer, Integer}) = Irrep(ir...)

Irrep(ir::Irrep) = ir

Base.show(io::IO, ir::Irrep) = print(io, "$(ir.l)$(ir.p == 1 ? 'e' : 'o')")

function Base.iterate(ir::Irrep, state::Int = 1)
    if state == 1
        return (ir.l, 2);
    elseif state == 2
        return (ir.p, 3);
    else
        return nothing;
    end
end

# Core collection traits
# Just a pragmatic approach
# Not sure if all of it is necessary
Base.length(::Irrep) = 2
Base.eltype(::Type{Irrep{T}}) where {T} = T
Base.first(ir::Irrep) = ir.l
Base.last(ir::Irrep) = ir.p

function Base.isless(a::Irrep, b::Irrep)
    key_a = (a.l, -a.p * (-1)^a.l)
    key_b = (b.l, -b.p * (-1)^b.l)
    return isless(key_a, key_b)
end

"""
    dim(x::Irreps)

Returns the dimension representation of `Irreps`, ``2l+1``
"""
dim(x::Irrep) = 2 * x.l + 1

"""
    isscalar(x::Irrep)

Check if the irrep is a scalar representation.
Equivalent to `x.l == 0 && x.p == 1`.
"""
isscalar(x::Irrep) = (x.l == 0) && (x.p == 1)

function iterator(::Type{Irrep}; lmax::Union{Int, Nothing} = nothing)
    return Channel{Irrep}() do ch
        for l in 0:typemax(Int)
            p1 = (-1)^l
            put!(ch, Irrep(l, p1))
            put!(ch, Irrep(l, -p1))
            if !isnothing(lmax) && l == lmax
                break
            end
        end
    end
end

function (T::Type{Irrep})(; lmax::Union{Int, Nothing} = nothing)
    iterator(T; lmax = lmax)
end

struct MulIrrep{T}
    mul::Int # might reconsider using T here, but not sure if it's worth the constraint
    ir::Irrep{T}
end

Base.show(io::IO, mul_ir::MulIrrep) = print(io, "$(mul_ir.mul)x$(mul_ir.ir)")

Base.isless(a::MulIrrep, b::MulIrrep) = isless((a.ir, a.mul), (b.ir, b.mul))

function Base.iterate(mul_ir::MulIrrep, state::Int = 1)
    if state == 1
        return (mul_ir.mul, 2);
    elseif state == 2
        return (mul_ir.ir, 3);
    else
        return nothing;
    end
end

dim(mul_ir::MulIrrep) = mul_ir.mul * dim(mul_ir.ir)

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

julia> Irrep("2e") in irs("0e + 2e")
true

# Empty Irreps
julia> Irreps(), Irreps("")
(, )
```
"""
struct Irreps{T}
    _irreps::Tuple{Vararg{MulIrrep{T}}}
end

Irreps(ir::Irrep{T}) where {T} = Irreps((MulIrrep(1, ir),))
Irreps(mul_ir::MulIrrep{T}) where {T} = Irreps((mul_ir,))
Irreps(mul_ir::Tuple{Int, Irrep{T}}) where {T} = Irreps((MulIrrep(mulir...),))
Irreps() = Irreps{Int}(tuple()) # Default to Int for empty constructor
Irreps(irs::Irreps{T}) where {T} = irs

function Irreps(irs::AbstractString)
    _irs = strip(irs)
    isempty(_irs) && return Irreps()
    mulirreps = Vector{MulIrrep{Int}}()
    if _irs != ""
        for mul_ir in split(_irs, "+")
            if occursin("x", mul_ir)
                mul, ir = split(mul_ir, "x")
                mul = parse(Int, mul)
                (mul >= 0) ||
                    throw(ArgumentError("mul should be greater than 0 got $mul"))
                irrep = Irrep(ir)
            else
                mul = 1
                irrep = Irrep(mul_ir)
            end
            push!(mulirreps, MulIrrep(mul, irrep))
        end
    end
    return Irreps(tuple(mulirreps...))
end

function Irreps(irs::Union{Vector, Tuple})
    isempty(irs) && return Irreps()

    # potentially mixed types
    initial_list = map(irs) do x
        if typeof(x) <: AbstractString
            # fix case for mixture of irrep and irrpes in a vector
            # if occursin("+", mul_irrep)
            #     irreps = Irreps(mul_irrep)
            # end
            MulIrrep(1, Irrep(x))
        elseif x isa MulIrrep
            x
        elseif x isa Tuple{Int, Irrep}
            MulIrrep(x...)
        elseif x isa Irrep
            MulIrrep(1, x)
        elseif length(x) == 2
            mul, ir = x
            MulIrrep(mul, Irrep(ir))
        else
            throw(ArgumentError("Invalid type in list for Irreps constructor: $(typeof(x))"))
        end
    end

    # common promoted type for the Irrep's T
    promoted_T = mapreduce(x -> eltype(x.ir), promote_type, initial_list)

    final_list = map(initial_list) do mul_ir
        new_ir = Irrep(
            convert(promoted_T, mul_ir.ir.l), convert(promoted_T, mul_ir.ir.p))
        MulIrrep(mul_ir.mul, new_ir)
    end

    return Irreps(tuple(final_list...))
end

function Base.show(io::IO, irs::Irreps)
    isempty(irs._irreps) && return print(io, "")
    print(io, join(string.(irs._irreps), "+"))
end

Base.iterate(irs::Irreps, args...) = iterate(irs._irreps, args...)
Base.getindex(irs::Irreps, i) = getindex(irs._irreps, i)
Base.length(irs::Irreps) = length(irs._irreps)
Base.eltype(irs::Irreps) = eltype(irs._irreps)
function Base.in(ir::Irrep, irs::Irreps)
    any(mul_ir -> mul_ir.ir == ir, irs._irreps)
end

function Base.:+(a::Irreps, b::Irreps)
    return Irreps((a._irreps..., b._irreps...))
end

function Base.:*(n::Integer, irs::Irreps)
    @assert n >= 0 "Multiplicity must be non-negative."
    return Irreps(tuple([MulIrrep(n * mul_ir.mul, mul_ir.ir) for mul_ir in irs]...))
end

Base.:*(irs::Irreps, n::Integer) = n * irs # Commutativity

function Base.:*(n::Integer, ir::Irrep)
    @assert n >= 0 "Multiplicity must be non-negative."
    return Irreps((MulIrrep(n, ir),))
end

Base.:+(ir1::Irrep, ir2::Irrep) = Irreps(ir1) + Irreps(ir2)

function Base.div(irs::Irreps, n::Integer)
    Irreps(tuple([MulIrrep(div(mul_ir.mul, mul_ir.ir), mul_ir.ir)
                  for mul_ir in irs]...))
end
function Base.fld(irs::Irreps, n::Integer)
    Irreps(tuple([MulIrrep(fld(mul_ir.mul, n), mul_ir.ir) for mul_ir in irs]...))
end

dim(irs::Irreps) = sum(dim, irs._irreps; init = 0)
num_irreps(irs::Irreps) = sum(x -> x.mul, irs._irreps; init = 0)
ls(irs::Irreps) = [mul_ir.ir.l for mul_ir in irs for _ in 1:mul_ir.mul]
function lmax(irs::Irreps)
    isempty(irs._irreps) ? -1 : maximum(x -> x.ir.l, irs._irreps)
end
function count(ir::Irrep, irs::Irreps)
    sum(mul_ir.mul for mul_ir in irs if mul_ir.ir == ir; init = 0)
end
function isscalar(irs::Irreps)
    all(mul_ir -> mul_ir.ir.l == 0 && mul_ir.ir.p == 1, irs._irreps)
end
function mul_gcd(irs::Irreps)
    isempty(irs._irreps) ? 0 : gcd(x -> x.mul, irs._irreps)
end

function Base.repeat(irs::Irreps, n::Integer)
    Irreps(tuple(repeat(collect(irs._irreps), n)...))
end

function remove_zero_multiplicities(irs::Irreps)
    Irreps(tuple([mul_ir for mul_ir in irs if mul_ir.mul > 0]...))
end

function unify(irs::Irreps)
    out = Vector{MulIrrep}{eltype(irs)}()
    for mul_ir in irs
        if !isempty(out) && last(out).ir == mul_ir.ir
            out[end] = MulIrrep(last(out).mul + mul_ir.mul, mul_ir.ir)
        else
            push!(out, mul_ir)
        end
    end
    return Irreps(tuple(out...))
end

simplify(irs::Irreps) = unify(remove_zero_multiplicities(irs))

function regroup(irs::Irreps)
    counts = Dict{Irrep, Int}()
    for mul_ir in irs
        counts[mul_ir.ir] = get(counts, mul_ir.ir, 0) + mul_ir.mul
    end
    sorted_irreps = sort(collect(keys(counts)))
    return Irreps(tuple([MulIrrep(counts[ir], ir)
                         for ir in sorted_irreps if counts[ir] > 0]...))
end

function Base.sort(irs::Irreps)
    v = collect(irs._irreps)
    p = sortperm(v)
    inv_p = invperm(p)
    return Irreps(tuple(v[p]...))
end

function Base.filter(f::Function, irs::Irreps)
    return Irreps(tuple([x for x in irs if f(x)]...))
end
