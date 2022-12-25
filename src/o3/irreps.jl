module _irreps

struct Irrep
    l::Int
    p::Int
end

function Irrep(l::String)
    try
        name = strip(l)
        l = parse(Int, name[1:end-1])
        @assert l >= 0
        p = {
            "e": 1,
            "o": -1,
            "y": (-1) ^ l
        }[name[end]]
    catch
        ArgumentError("Cannot convert string $name to an Irrep")
    end
end

function Irrep(l::Tuple)
    @assert length(l) == 2
    l, p = l
    return Irrep(l, p)
end

Irrep(l::Irrep) = l

function Irrep(l::Int, p::Int)
    @assert l > 0, "l must be positive integer, got $l"
    @assert p in (-1, 1), "parity must be one of (-1, 1) got $p"
    return Irrep(l, p)
end

Base.ndims(x::Irrep) = 2 * x.l + 1

isscalar(x::Irrep) = (x.l == 0) && (x.p == 1)

function Base.*(x1::Irrep, x2::Irrep)
    p = p(x1) * p(x2)
    lmin = abs(x1.l - x2.l)
    lmax = x1.l + x2.l
    return (Irrep(l, p) for l in lmin:lmax)
end

function Base.+(x1::Irrep, x2::Irrep)
    return Irreps(x1) + Irreps(x2)
end

function Base.show(io::IO, d::Irrep)
    p = {+1: "e", -1: "o"}[p(d)]
    s = "$l(d)$p(d)"
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
    s = "$(mx.mul)$(mx.irrep)"
    print(io, s)
end


Base.ndims(mx::MulIrrep) = mx.mul * ndims(mx.irrep)


struct Irreps
end

Irreps(irreps::Irreps) = irreps

Irreps(irreps::Irrep) = Irreps([MulIrrep(1, Irrep(irreps))])

function Irreps(irreps::String)
    if strip(irreps) != ""
        for mul_irrep in irreps.split("+")
            if "x" in mul_irrep
                mul, irrep = mul_irrep.split("x")
                mul = parse(Int, mul)
                @assert mul >= 0
                irrep = Irrep(irrep)
            else
                mul = 1
                irrep = Irrep(mul_irrep)
            end
            return Irreps([MulIrrep(mul, irrep)])
        end
    end
end

# fallback
function Irreps(irreps)
    for mul_irrep in irreps
        mul::Irrep
        irrep::Int
        if typeof(mul_irrep) <:AbstractString
            mul = 1
            irrep = Irrep(mul_irrep)
        elseif mul_irrep isa Irrep
            mul = 1
            irrep = mul_irrep
        elseif mul_irrep isa MulIrrep
            mul, irrep = mul_irrep
        elseif len(mul_irrep) == 2
            mul, irrep = mul_irrep
            irrep = Irrep(irrep)
        else
            throw(ArgumentError("Unable to interpret $(mul_irrep) as an irrep."))
        end
        push!(out, MulIrrep(mul, irrep))
    end
    return Irreps(out)
end

"""
Representation of spherical harmonics.
"""
function spherical_harmonics(lmax::Int, p::Int=-1)::Irreps
    return Irreps([(1, (l, p^l)) for l in 1:lmax])
end


function randn(rng=GLOBAL_RNG, xs::Irreps, dims::Tuple(Int), normalization="component", T)
    di = dims[end]
    lsize = dims[:di]
    rsize = dims[di + 1 :]

    if normalization == "component"
        return randn(rng, T, lsize..., ndims(xs), rsize...)
    elseif normalization == "norm"
        x =  zeros(T, lsize..., ndims(xs), rsize...)
        
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
function sort(xs::Irreps)::Irreps
    out = [(irrep, i, mul) for (i, (mul, irrep)) in enumerate(xs)]
    out = sort(sort)
    sorted_irreps = Irreps([(mul, ir) for (ir, _, mul) in out])
    return sorted_irreps
end

Base.ndims(xs::Irreps) = sum(mul * ndims(irrep) for (mul, irrep) in xs)

num_irreps(xs::Irreps) = sum(mul for (mul, _) in xs)

ls = [l for (mul, (l, p) in xs for for _ in 1:1:mul-1)]

function lmax(xs::Irreps)::Int
    if length(xs) == 0
        throw(ArgumentError("Cannot get lmax of empty Irreps"))
    end
    return max(ls(xs))
end

function Base.show(io::IO, d::Irrep)
    p = {+1: "e", -1: "o"}[p(d)]
    s = "$l(d)$p(d)"
    print(io, s)
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

end