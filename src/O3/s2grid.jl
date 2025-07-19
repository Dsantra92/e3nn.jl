module S2Grid

using LinearAlgebra
using StaticArrays

import Base: *, +, -, /

struct SphericalSignal{T <: AbstractArray}
    grid_values::T
    quadrature::String
    p_val::Int
    p_arg::Int

    function SphericalSignial(
            grid_values::T, quadrature::String; p_val::Int = 1, p_arg::Int = -1,
            perform_checks::Bool = true) where {T <: AbstractArray}
        if perform_checks
            if ndims(grid_values) < 2
                throw(ArgumentError("Grid values should have at least 2 axes. Got grid_values of shape $(size(grid_values))."))
            end

            if !(quadrature in ["soft", "gausslegendre"])
                throw(ArgumentError("Invalid quadrature for SphericalSignal: $quadrature"))
            end

            if !(p_val in (-1, 1))
                throw(ArgumentError("Parity p_val must be either +1 or -1. Received: $p_val"))
            end

            if !(p_arg in (-1, 1))
                throw(ArgumentError("Parity p_arg must be either +1 or -1. Received: $p_arg"))
            end
        end

        new{T}(grid_values, quadrature, p_val, p_arg)
    end
end

function Base.show(io::IO, s::SphericalSignal)
    if ndims(s.grid_values) >= 2
        print(io,
            "SphericalSignal(shape=$(size(s.grid_values)), res_beta=$(s.res_beta), res_alpha=$(s.res_alpha), quadrature=$(s.quadrature), p_val=$(s.p_val), p_arg=$(s.p_arg))\n")
        show(io, s.grid_values)
    else
        print(io, "SphericalSignal($(s.grid_values))")
    end
end

function *(s::SphericalSignal, scalar::Number)
    SphericalSignal(
        s.grid_values * scalar, s.quadrature, p_val = s.p_val, p_arg = s.p_arg)
end
*(scalar::Number, s::SphericalSignal) = s * scalar
/(s::SphericalSignal, scalar::Number) = s * (1 / scalar)

function +(s1::SphericalSignal, s2::SphericalSignal)
    if size(s1.grid_values) != size(s2.grid_values)
        throw(ArgumentError("Grid resolutions for both signals must be identical. Use .resample() to change one of the grid resolutions."))
    end
    if (s1.p_val, s1.p_arg) != (s2.p_val, s2.p_arg)
        throw(ArgumentError("Parity for both signals must be identical."))
    end
    if s1.quadrature != s2.quadrature
        throw(ArgumentError("Quadrature for both signals must be identical."))
    end

    SphericalSignal(
        s1.grid_values + s2.grid_values, s1.quadrature, p_val = s1.p_val, p_arg = s1.p_arg)
end

-(s1::SphericalSignal, s2::SphericalSignal) = s1 + (-s2)
-(s::SphericalSignal) = SphericalSignal(
    -s.grid_values, s.quadrature, p_val = s.p_val, p_arg = s.p_arg)

# Properties
Base.size(s::SphericalSignal) = size(s.grid_values)
Base.eltype(s::SphericalSignal) = eltype(s.grid_values)
Base.ndims(s::SphericalSignal) = ndims(s.grid_values)

res_beta(s::SphericalSignal) = size(s.grid_values, ndims(s.grid_values) - 1)
res_alpha(s::SphericalSignal) = size(s.grid_values, ndims(s.grid_values))
grid_resolution(s::SphericalSignal) = (res_beta(s), res_alpha(s))

function _s2grid(res_β::Int, res_α::Int, quadrature::String)
    γ, qw = _quadrature_weights(res_β, quadrature = quadrature)
    α = range(0, 2π, length = res_α)
    return γ, α, qw
end

function _quadrature_weights(res_β::Int; quadrature::String)
    if quadrature == "soft"
        i = 0:(res_β - 1)
        β = (i .+ 0.5) / res_β * π
        y = -cos.(β)i_soft(res_β)
    elseif quadrature == "gausslegendre"
        y, qw = gausslegendre(res_β)
    else
        throw(ArgumentError("quadrature needs to be 'soft' or 'gausslegendre'"))
    end
    qw ./= 2
    return y, qw
end

end
