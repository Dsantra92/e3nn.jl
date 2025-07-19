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

Instantiate a new [`o3.Irrep`](@ref) object from a tuple of (`l`, `p`).
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
