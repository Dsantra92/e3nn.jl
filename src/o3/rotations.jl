module rot

using LinearAlgebra
using Quaternions
using Rotations
using StaticArrays

export Quaternion
export euler_angles, CartesianToSphericalAngles, SphercialAnglesToCartesian

# Cross conversion

function (::Type{Q})(aa::AngleAxis) where {Q <: Quaternion}
    s, c = sincos(aa.theta / 2)
    return Q(c, s * aa.axis_x, s * aa.axis_y, s * aa.axis_z)
end

function (::Type{AA})(q::Quaternion) where {AA <: AngleAxis}
    AA(QuatRotation(q, false)) # already implemented method
end

function (::Type{AA})(; α::Real, β::Real, γ::Real) where {AA <: AngleAxis}
    AA(RotYXY(promote(α, β, γ)...))
end

function (::Type{Q})(; α::Real, β::Real, γ::Real) where {Q <: Quaternion}
    return QuatRotation(RotYXY(promote(α, β, γ)...)).q |> Q
end

function euler_angles(x::T) where {T <: Union{QuatRotation, RotMatrix3, AngleAxis}}
    Rotations.params(RotYXY(x))
end

euler_angles(q::Q) where {Q <: Quaternion} = euler_angles(QuatRotation(q, false))

function (::Type{R})(q::Quaternion) where {R <: RotMatrix3}
    QuatRotation(q) |> R
end

function (::Type{Q})(R::RotMatrix3) where {Q <: Quaternion}
    return QuatRotation(R).q |> Q
end

function CartesianToSphericalAngles(x::SVector{3, T}) where {T <: Real}

    # done in e3nn to remove NaNs
    # need to check for Julia
    normalize!(x, 2)
    clamp!(x, -1, 1)
    return SVector(acos[2], atan(x[1], x[3]))
end

function SphercialAnglesToCartesian(α::Real, β::Real)
    α, β = promote(α, β)
    sβ, cβ = sincos(β)
    sα, cα = sincos(α)
    return SVector(sβ * sα, cβ, sβ * cα)
end

end
