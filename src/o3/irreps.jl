using Random: GLOBAL_RNG, AbstractRNG

struct Irrep
    l::Int
    p::Int

    function Irrep(l::Int, p::Int)
        (l >= 0) || throw(ArgumentError("l must be zero or positive integer, got $l"))
        (p in (-1, 1)) || throw(ArgumentError("parity(p) must be one of (-1, 1) got $p"))
        return new(l, p)
    end
end

function Irrep(l::T) where T<:AbstractString
    name = strip(l)
    try
        l = parse(Int, name[1:end-1])
        (l >= 0) || throw(ArgumentError("l must be zero or positive integer, got $l"))
        p = Dict(
            'e' => 1,
            'o' => -1,
            'y' => (-1) ^ l
        )[name[end]]
        return Irrep(l, p)
    catch
        throw(ArgumentError("Cannot convert string $name to an Irrep"))
    end
end

function Irrep(l::Tuple)
    @assert length(l) == 2
    l, p = l
    return Irrep(l, p)
end

Irrep(l::Irrep) = l

dim(x::Irrep) = 2 * x.l + 1

isscalar(x::Irrep) = (x.l == 0) && (x.p == 1)

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

# Used for comparison
Base.isless(x1::Irrep, x2::Irrep) = Base.isless((x1.l, x1.p), (x2.l, x2.p)) 

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

dim(mx::MulIrrep) = mx.mul * dim(mx.irrep)

Base.convert(::Type{Tuple}, x::MulIrrep) = (x.mul, x.irrep)
struct Irreps
    irreps::Vector{MulIrrep}
end

Irreps(irrep::Irrep) = Irreps([MulIrrep(1, irrep)])
Irreps(irreps::Irreps) = irreps

function Irreps(irreps::T) where T<:AbstractString
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
        if typeof(mul_irrep) <:AbstractString
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

Base.in(x::Irrep, xs::Irreps) = x ∈ [mx.irrep for mx in xs.irreps]
Base.count(x::Irrep, xs::Irreps) = sum([mx.mul for mx in xs.irreps if mx.irrep == x], init=0)
Base.iterate(xs::Irreps, state=1) = state > length(xs) ? nothing : (xs[state], state+1)

"""
Representation of spherical harmonics.
"""
function spherical_harmonics(lmax::Int, p::Int=-1)::Irreps
    return Irreps([(1, (l, p^l)) for l in 1:lmax])
end

function Base.randn(xs::Irreps, dims, normalization::String, rng::AbstractRNG, ::Type{T}) where T<:Number
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
    for (out, irreps) in xs
        if out && out[end][2] == Irrep
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
remove_zero_multiplicities(xs::Irreps) = [(mul, irreps) for (mul, irreps) in xs if mul > 0] |> Irreps

"""
Sort the representations.
"""
function Base.sort(xs::Irreps)::Irreps
    out = [(mx.irrep, i, mx.mul) for (i, mx) in enumerate(xs)]
    out = Base.sort(out)
    sorted_irreps = Irreps([(mul, irrep) for (irrep, _, mul) in out])
    return sorted_irreps
end

dim(xs::Irreps) = sum([mx.mul * dim(mx.irrep) for mx in xs], init=0)

num_irreps(xs::Irreps) = sum([mx.mul for mx in xs], init=0)

ls(xs::Irreps) = [mx.irrep.l for mx in xs for _ in 1:mx.mul]

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
