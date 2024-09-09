using LinearAlgebra
using StaticArrays

function su2_generators(j::Int)
    m = range(-j, j - 1, step = 1)
    raising = diagm(-1 => -sqrt.(j * (j + 1) .- m .* (m .+ 1)))

    m = range(-j + 1, j, step = 1)
    lowering = diagm(1 => sqrt.(j * (j + 1) .- m .* (m .- 1)))

    m = range(-j, j, step = 1)
    return stack(
        [0.5 * (raising + lowering),             # x (usually)
            Diagonal(1im * m),                   # z (usually)
            -0.5im * (raising - lowering)],      # -y (usually)
        dims = 3)
end

function change_basis_real_to_complex(l::Int)
    # https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
    q = zeros(ComplexF64, 2 * l + 1, 2 * l + 1)
    for m in (-l):-1
        q[l + m + 1, l + abs(m) + 1] = 1 / sqrt(2)
        q[l + m + 1, l - abs(m) + 1] = -im / sqrt(2)
    end
    q[l + 1, l + 1] = 1
    for m in 1:l
        q[l + m + 1, l + abs(m) + 1] = (-1)^m / sqrt(2)
        q[l + m + 1, l - abs(m) + 1] = im * (-1)^m / sqrt(2)
    end
    q = (-im)^l * q  # Added factor of im^l to make the Clebsch-Gordan coefficients real
    return q
end

function so3_generators(l::Int)
    X = su2_generators(l)
    Q = change_basis_real_to_complex(l)
    Q_c_T = conj(transpose(Q))

    for i in 1:size(X, 3)
        @views X[:, :, i] = Q_c_T * X[:, :, i] * Q
    end
    return real.(X)
end

function wigner_D(l::Int, α::T, β::T, γ::T) where {T <: Real}
    X = so3_generators(l)
    return exp(α * X[:, :, 2]) * exp(β * X[:, :, 1]) * exp(γ * X[:, :, 2])
end

function wigner_D(l::Int, angles::Tuple{Vararg{T, 3}}) where {T <: Real}
    return wigner_D(l, angles[1], angles[2], angles[3])
end

function wigner_D(l::Int, angles::SVector{3, T}) where {T <: Real}
    return wigner_D(l, angles[1], angles[2], angles[3])
end

Broadcast.broadcast(::typeof(wigner_D), l, α, β, γ) = broadcast(wigner_D, Ref(l), α, β, γ)
