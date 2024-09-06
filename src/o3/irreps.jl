using Random: AbstractRNG
"""
Irrreducible representations of ``O(3)``.

This struct does not contain any data; it is a structure that describes the representation.
It is typically used as an argument in other parts of the library to define the input and output representations of
functions.

# Fields:
- `l::Int`: non-negative integer, the degree of the representation, ``l = 0, 1, \\dots``
- `p::Int`: the parity of the representation, either 1 (even) or -1 (odd)
"""
struct Irrep
    l::Int
    p::Int

    @doc """
         Irrep(l::Int, p::Int)

     Instantiate a new [`o3.Irrep`](@ref) object.

     # Arguments:
     - `l::Int`: non-negative integer, the degree of the representation, ``l = 0, 1, \\dots``
     - `p::Int`: the parity of the representation, either ``1`` (even) or ``-1`` (odd)

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
    function Irrep(l::Int, p::Int)
        (l >= 0) || throw(ArgumentError("l must be zero or positive integer, got $l"))
        (p in (-1, 1)) || throw(ArgumentError("parity(p) must be one of (-1, 1) got $p"))
        return new(l, p)
    end
end

"""
    Irrep(ir::T) where {T <: AbstractString}

Instantiate a new [`o3.Irrep`](@ref) object from it's string representation.

# Arguments
- `s::String`: A string representation of the irrep.

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
function Irrep(ir::T) where {T <: AbstractString}
    name = strip(ir)
    try
        l = parse(Int, name[1:(end - 1)])
        (l >= 0) || throw(ArgumentError("l must be zero or positive integer, got $l"))
        p = Dict('e' => 1, 'o' => -1, 'y' => (-1)^l)[name[end]]
        return Irrep(l, p)
    catch
        throw(ArgumentError("Cannot convert string $name to an Irrep"))
    end
end

"""
    Irrep(ir::Tuple{Int, Int})

Instantiate a new [`o3.Irrep`](@ref) object from a tuple of (`l`, `p`).
"""
Irrep(ir::Tuple{Int, Int}) = Irrep(ir...)

Irrep(l::Irrep) = l

dim(x::Irrep) = 2 * x.l + 1

"""
    isscalar(x::Irrep)

Check if the irrep is a scalar representation.
Equivalent to `x.l == 0 && x.p == 1`.
"""
isscalar(x::Irrep) = (x.l == 0) && (x.p == 1)

Base.iterate(x::Irrep, args...) = iterate((x.l, x.p), args...)

function Base.:*(x1::Irrep, x2::Irrep)
    p = x1.p * x2.p
    lmin = abs(x1.l - x2.l)
    lmax = x1.l + x2.l
    return (Irrep(l, p) for l in lmin:lmax)
end

function Base.:*(i::Int, x::Irrep)
    return Irreps([(i, x)])
end

function Base.:+(x1::Irrep, x2::Irrep)
    return Irreps(x1) + Irreps(x2)
end

function Base.isless(x1::Irrep, x2::Irrep)
    Base.isless((x1.l, -x1.p * (-1)^x1.l), (x2.l, -x2.p * (-1)^x2.l))
end

function Base.show(io::IO, x::Irrep)
    p = Dict(+1 => "e", -1 => "o")[x.p]
    s = "$(x.l)$p"
    print(io, s)
end

function D_from_angles(x::Irrep, α, β, γ, k)
    throw(error("Not Implemnted yet"))
end

function D_from_quaternion(x::Irrep, q, k)
    throw(error("Not Implemnted yet"))
end

function D_from_matrix(x::Irrep, R)
    throw(error("Not Implemnted yet"))
end

function D_from_matrix(x::Irrep, axis, angle)
    throw(error("Not Implemnted yet"))
end

struct MulIrrep
    mul::Int
    irrep::Irrep
end

function Base.show(io::IO, mx::MulIrrep)
    s = "$(mx.mul)x$(mx.irrep)"
    print(io, s)
end

"""
dim(x::MulIrrep)
"""
dim(mx::MulIrrep) = mx.mul * dim(mx.irrep)

Base.convert(::Type{Tuple}, x::MulIrrep) = (x.mul, x.irrep)

struct Irreps
    irreps::Vector{MulIrrep}
end

Irreps(irrep::Irrep) = Irreps([MulIrrep(1, irrep)])
Irreps(irreps::Irreps) = irreps

function Irreps(irreps::T) where {T <: AbstractString}
    mulirreps = MulIrrep[]
    if strip(irreps) != ""
        for mul_irrep in split(irreps, "+")
            if occursin("x", mul_irrep)
                mul, irrep = split(mul_irrep, "x")
                mul = parse(Int, mul)
                (mul >= 0) || throw(ArgumentError("mul should be greater than 0 got $mul"))
                irrep = Irrep(irrep)
            else
                mul = 1
                irrep = Irrep(mul_irrep)
            end
            push!(mulirreps, MulIrrep(mul, irrep))
        end
    end
    return Irreps(mulirreps)
end

# fallback
function Irreps(irreps::AbstractVector)
    mul = nothing
    irrep = nothing
    out = MulIrrep[]
    for mul_irrep in irreps
        if typeof(mul_irrep) <: AbstractString
            # fix case for mixture of irrep and irrpes in a vector
            # if occursin("+", mul_irrep)
            #     irreps = Irreps(mul_irrep)
            # end
            mul = 1
            irrep = Irrep(mul_irrep)
        elseif mul_irrep isa Irrep
            mul = 1
            irrep = mul_irrep
        elseif mul_irrep isa MulIrrep
            mul, irrep = mul_irrep
        elseif length(mul_irrep) == 2
            mul, irrep = mul_irrep
            irrep = Irrep(irrep)
        else
            throw(ArgumentError("Unable to interpret $(mul_irrep) as an irrep."))
        end
        push!(out, MulIrrep(mul, irrep))
    end
    return Irreps(out)
end

function Base.show(io::IO, xs::Irreps)
    s = join(["$(mul_irrep)" for mul_irrep in xs.irreps], "+")
    print(io, s)
end

Base.:(==)(x1::Irreps, x2::Irreps) = x1.irreps == x2.irreps
Base.:+(x1::Irreps, x2::Irreps) = Irreps(vcat(x1.irreps, x2.irreps))
Base.:*(i::Int, xs::Irreps) = Irreps(repeat(xs.irreps, i))
Base.:*(xs::Irreps, i::Int) = Base.:*(i::Int, xs::Irreps)

Base.length(xs::Irreps) = length(xs.irreps)

Base.getindex(xs::Irreps, i::Int) = xs.irreps[i]
Base.getindex(xs::Irreps, idx::AbstractRange) = xs.irreps[idx] |> Irreps
Base.firstindex(xs::Irreps) = Base.firstindex(xs.irreps)
Base.lastindex(xs::Irreps) = Base.lastindex(xs.irreps)

# TODO: implement indexing by mul and dim if required

Base.in(x::Irrep, xs::Irreps) = x ∈ [mx.irrep for mx in xs.irreps]
function Base.count(x::Irrep, xs::Irreps)
    sum([mx.mul for mx in xs.irreps if mx.irrep == x], init = 0)
end
function Base.iterate(xs::Irreps, state = 1)
    state > length(xs) ? nothing : (xs[state], state + 1)
end

"""
Representation of spherical harmonics.
"""
function spherical_harmonics(lmax::Int, p::Int = -1)::Irreps
    return Irreps([(1, (l, p^l)) for l in 1:lmax])
end

function Base.randn(
        xs::Irreps,
        dims,
        normalization::String,
        rng::AbstractRNG,
        ::Type{T}
) where {T <: Number}
    di = dims[end]
    lsize = dims[:di]
    rsize = dims[di + 1 :]

    if normalization == "component"
        return randn(rng, T, lsize..., dim(xs), rsize...)
        # implement the norm
    end
end

"""
Simplify the representaions.
"""
function simplify(xs::Irreps)::Irreps
    out = []
    xs = Base.sort(xs)
    for mul_ir in xs
        mul, irrep = mul_ir.mul, mul_ir.irrep
        if length(out) != 0 && out[end][2] == irrep
            out[end] = (out[end][1] + mul, irrep)
        elseif mul > 0
            push!(out, (mul, irrep))
        end
    end
    return Irreps(out)
end

"""
Remove any irreps with multiplicities of zero.
"""
remove_zero_multiplicities(xs::Irreps) = [(mul, irreps) for (mul, irreps) in xs if mul > 0] |>
                                         Irreps

"""
    Base.sort(xs::Irreps)

Sort the representations.

# Examples
```julia-repl
julia> sort(Irreps([(1, (1, 1)), (1, (2, 1)), (1, (1, -1))]))
```
"""
function Base.sort(xs::Irreps)
    out = [(mx.irrep, i, mx.mul) for (i, mx) in enumerate(xs)]
    out = Base.sort(out)
    sorted_irreps = Irreps([(mul, irrep) for (irrep, _, mul) in out])
    return sorted_irreps
end

dim(xs::Irreps) = sum([mx.mul * dim(mx.irrep) for mx in xs], init = 0)

num_irreps(xs::Irreps) = sum([mx.mul for mx in xs], init = 0)

ls(xs::Irreps) = [mx.irrep.l for mx in xs for _ in 1:(mx.mul)]

"""
    Base.iterate(::Type{Irrep}, state=(0, true))

Iterator through all the irreps of O(3).

# Examples
```julia-repl
julia> first(Irrep)
0e

julia> collect(Iterators.take(Irrep, 6)) # set lmax as 6
6-element Vector{Any}:
 0e
 0o
 1o
 1e
 2e
 2o
```
"""
function Base.iterate(::Type{Irrep}, state = (0, true))
    l, is_positive = state
    if is_positive
        return (Irrep(l, (-1)^l), (l, false))
    else
        return (Irrep(l, -(-1)^l), (l + 1, true))
    end
end

Base.IteratorSize(::Type{Irrep}) = Base.IsInfinite()
Base.IteratorEltype(::Type{Irrep}) = Base.HasEltype()
Base.eltype(::Type{Irrep}) = Irrep

function lmax(xs::Irreps)::Int
    if length(xs) == 0
        throw(ArgumentError("Cannot get lmax of empty Irreps"))
    end
    return maximum(ls(xs))
end

function D_from_angles(xs::Irreps, α, β, γ, k)
    throw(error("Not Implemnted yet"))
end

function D_from_quaternion(xs::Irreps, q, k)
    throw(error("Not Implemnted yet"))
end

function D_from_matrix(xs::Irreps, R)
    throw(error("Not Implemnted yet"))
end

function D_from_axis_angle(xs::Irreps, axis, angle)
    throw(error("Not Implemnted yet"))
end
