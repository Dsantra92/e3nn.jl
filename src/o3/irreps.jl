using Random: AbstractRNG
using Rotations
using Rotations: params
using Quaternions
using BlockDiagonals: BlockDiagonal

"""
using StaticArrays: AngleAxis
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
         Irrep(l::Integer, p::Integer)

     Instantiate a new [`o3.Irrep`](@ref) object.

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
    Irrep(ir::Tuple{Integer, Integer})

Instantiate a new [`o3.Irrep`](@ref) object from a tuple of (`l`, `p`).
"""
Irrep(ir::Tuple{Integer, Integer}) = Irrep(ir...)

Irrep(l::Irrep) = l

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

Base.iterate(x::Irrep, args...) = iterate((x.l, x.p), args...)

"""
    Base.:*(x1::Irrep, x2::Irrep)

Returns a generator from the product of 2 `Irrep`.

```jldoctest
julia> generator = Irrep("1o") * Irrep("2y");

julia> collect(generator)
3-element Vector{Irrep}:
 1o
 2o
 3o

```
"""
function Base.:*(x1::Irrep, x2::Irrep)
    p = x1.p * x2.p
    lmin = abs(x1.l - x2.l)
    lmax = x1.l + x2.l
    return (Irrep(l, p) for l in lmin:lmax)
end

function Base.:*(i::Integer, x::Irrep)
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

"""
    wigner_D(x::Irrep, α::T, β::T, γ::T, k::T = 0.0) where {T <: Real}

Returns the WingerD matrix for the given `Irrep`, euler angles(α, β, γ) and k.
"""
function wigner_D(x::Irrep, α::T, β::T, γ::T, k::T = 0.0) where {T <: Real}
    return wigner_D(x.l, α, β, γ) .* x.p^k
end

"""
    wigner_D(x::Irrep, q::Quaternion, k::T) where {T <: Real}

Returns the WingerD matrix for the given `Irrep`, `Quaternion` and k.
"""
function wigner_D(x::Irrep, q::Quaternion, k::T) where {T <: Real}
    α, β, γ = q |> QuatRotation |> RotYXY |> params
    return wigner_D(x, α, β, γ, k)
end

"""
    wigner_D(x::Irrep, q::Quaternion, k::T) where {T <: Real}

Returns the WingerD matrix for the given `Irrep`, `Rotations.RotMatrix3` and k.
"""
function wigner_D(x::Irrep, R::RotMatrix3{T}) where {T <: Real}
    R = RotYXY(R)
    d = R |> det |> sign
    R .*= d
    k = (1 - d) / 2
    return wigner_D(x, params(R)..., k)
end

"""
    wigner_D(x::Irrep, q::Quaternion, k::T) where {T <: Real}

Returns the WingerD matrix for the given `Irrep`, `Rotations.AngleAxis` and k.
"""
function wigner_D(x::Irrep, aa::AngleAxis)
    R = RotYXY(aa)
    return wigner_D(x, R)
end

struct MulIrrep
    mul::Int
    irrep::Irrep
end

MulIrrep(mul::Integer, irrep::AbstractString) = MulIrrep(mul, Irrep(irrep))

function Base.show(io::IO, mx::MulIrrep)
    s = "$(mx.mul)x$(mx.irrep)"
    print(io, s)
end

dim(mx::MulIrrep) = mx.mul * dim(mx.irrep)

Base.:(==)(x1::MulIrrep, x2::Tuple) = (x1.mul == x2[1]) && (x1.irrep == x2[2])

"""
    Irreps

Direct sum of irreducible representations of ``O(3)``.

This struct does not contain any data, it is a structure that describes the representation.
It is typically used as an argument of other types in the library to define the input and output representations of
functions.

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
    irreps::Vector{MulIrrep}
end

Irreps(irrep::Irrep) = Irreps([MulIrrep(1, irrep)])
Irreps(irreps::Irreps) = irreps
Irreps() = Irreps([])

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
Base.:*(i::Integer, xs::Irreps) = Irreps(repeat(xs.irreps, i))
Base.:*(xs::Irreps, i::Integer) = Base.:*(i::Integer, xs::Irreps)

Base.length(xs::Irreps) = length(xs.irreps)

Base.getindex(xs::Irreps, i::Integer) = xs.irreps[i]
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
    multiplicity(xs::Irreps, ir::Irrep)

Multiplicity of the given `Irrep` in the `Irreps`.
"""
multiplicity(xs::Irreps, ir::Irrep) = [mul for (; mul, irrep) in xs if irrep == ir] |> sum

"""
    spherical_harmonics(lmax::Integer, p::Integer = -1)

Representation of spherical harmonics for given `lmax` and parity.

# Examples
```jldocket
julia> spherical_harmonics(3)
1x0e+1x1o+1x2e+1x3o

julia> spherical_harmonics(4, p=1)
1x0e+1x1e+1x2e+1x3e+1x4e
```
"""
function spherical_harmonics(lmax::Integer, p::Integer = -1)
    return Irreps([(1, (l, p^l)) for l in 1:lmax])
end

"""
    simplify(xs::Irreps)

Simplify the representaions.
Simplify does not sort the representations, so equivalent representations which
are separated from each other are not combined.
Use `regroup`(@ref) instead for such use cases.

# Examples
julia> Irreps("1e + 1e + 0e") |> simplify
2x1e+1x0e

julia> Irreps("1e + 1e + 0e + 1e") |> simplify # no sorting
2x1e+1x0e+1x1e
"""
function simplify(xs::Irreps)::Irreps
    out = []
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
    unify(self::Irreps)

Regroup same irreps together.

# Examples
```jldoctest
julia> unify(Irreps("0e + 1e"))
1x0e+1x1e

julia> xs = Irreps("0e + 1e + 1e")

julia> unify(xs)
1x0e+2x1e

julia> unify(Irreps("0e + 0x1e + 0e"))
1x0e+0x1e+1x0e
```
"""
function unify(xs::Irreps)
    out = Tuple{Int, Irrep}[]
    for (; mul, irrep) in xs
        if !isempty(out) && out[end][2] == irrep
            out[end] = (out[end][1] + mul, irrep)
        else
            push!(out, (mul, irrep))
        end
    end
    return Irreps(out)
end

"""
    regroup(xs::Irreps)

Regroup the same irrep together.
Equivalent to `sort` followed by `simplify`.

# Examples
julia> Irreps("1e + 0e + 1e + 0x2e") |> regroup
1x0e+2x1e

"""
regroup(xs::Irreps) = xs |> sort |> simplify

"""
    remove_zero_multiplicities(xs::Irreps)

Remove any irreps with multiplicities of zero.
# Examples
```jldoctest
julia> Irreps("4x0e + 0x1o + 2x3e") |> remove_zero_multiplicities
4x0e+2x3e
```
"""
remove_zero_multiplicities(xs::Irreps) = [(mul, irrep)
                                          for (; mul, irrep) in xs if mul > 0] |>
                                         Irreps

"""
    Base.sort(xs::Irreps)

Sort the representations.

# Examples
```jldoctest
julia> sort(Irreps([(1, (1, 1)), (1, (2, 1)), (1, (1, -1))]))
1x1o+1x1e+1x2e
```
"""
function Base.sort(xs::Irreps)
    out = [(mx.irrep, i, mx.mul) for (i, mx) in enumerate(xs)]
    out = Base.sort(out)
    sorted_irreps = Irreps([(mul, irrep) for (irrep, _, mul) in out])
    return sorted_irreps
end

"""
    dim(xs::Irreps)

Dimension of the irreps.
# Examples:
```
julia> Irreps("3x0e + 2x1e") |> dim
9
```
"""
dim(xs::Irreps) = sum([mx.mul * dim(mx.irrep) for mx in xs], init = 0)

"""
    num_irreps(xs::Irreps)

Sum of the multiplicities.

# Examples:
```
julia> Irreps("3x0e + 2x1e") |> num_irreps
5

```
"""
num_irreps(xs::Irreps) = sum([mx.mul for mx in xs], init = 0)

"""
    ls(xs::Irreps)

List of the l values.

# Examples:
```jldoctest
julia> Irreps("3x0e + 2x1e") |> ls
5-element Vector{Int64}:
 0
 0
 0
 1
 1

```
"""
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

"""
    lmax(xs::Irreps)

Returns the maximum value of l in the Irreps.

# Examples
```jldoctest
julia> Irreps("3x0e + 2x1e") |> lmax
1
```
"""
function lmax(xs::Irreps)
    if length(xs) == 0
        throw(ArgumentError("Cannot get lmax of empty Irreps"))
    end
    return maximum(ls(xs))
end

function wigner_D(xs::Irreps, α::T, β::T, γ::T, k::T = 0.0) where {T <: Real}
    return BlockDiagonal([wigner_D(irrep, α, β, γ, k) for (; mul, irrep) in xs
                          for _ in 1:mul])
end

function wigner_D(xs::Irreps, q::Quaternion, k::T) where {T <: Real}
    α, β, γ = q |> QuatRotation |> RotYXY |> params
    return wigner_D(xs, α, β, γ, k)
end

function wigner_D(xs::Irreps, R::RotMatrix3{T}) where {T <: Real}
    R = RotYXY(R)
    d = R |> det |> sign
    R .*= d
    k = (1 - d) / 2
    return wigner_D(xs, params(R)..., k)
end

function wigner_D(xs::Irreps, aa::AngleAxis)
    R = RotYXY(aa)
    return wigner_D(xs, R)
end
