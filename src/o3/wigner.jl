using LinearAlgebra
using StaticArrays
using Quaternions
using TensorOperations

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

function su2_clebsch_gordan(j1::T, j2::T, j3::T) where {T <: Real}
    mat = zeros(Float64, Int(2 * j1 + 1), Int(2 * j2 + 1), Int(2 * j3 + 1))

    if Int(2 * j3) in Int(2 * abs(j1 - j2)):2:Int(2 * (j1 + j2))
        for m1 in [x / 2 for x in (-Int(2 * j1)):2:Int(2 * j1)]
            for m2 in [x / 2 for x in (-Int(2 * j2)):2:Int(2 * j2)]
                if abs(m1 + m2) <= j3
                    val = _su2_cg(promote(j1, m1), promote(j2, m2), promote(j3, m1 + m2))
                    mat[Int(j1 + m1 + 1), Int(j2 + m2 + 1), Int(j3 + m1 + m2 + 1)] = val
                end
            end
        end
    end
    return mat / sqrt(2 * j3 + 1)
end

function f(n::T) where {T <: Real}
    return factorial(round(Int, n))
end

function _su2_cg(idx1::Tuple{Float64, Float64}, idx2::Tuple{Float64, Float64},
        idx3::Tuple{Float64, Float64})
    j1, m1 = idx1
    j2, m2 = idx2
    j3, m3 = idx3

    if m3 != m1 + m2
        return 0.0
    end

    vmin = Int(max(-j1 + j2 + m3, -j1 + m1, 0))
    vmax = Int(min(j2 + j3 + m1, j3 - j1 + j2, j3 + m3))

    C = sqrt((2.0 * j3 + 1.0) *
             (f(j3 + j1 - j2) * f(j3 - j1 + j2) *
              f(j1 + j2 - j3) * f(j3 + m3) *
              f(j3 - m3)) //
             (f(j1 + j2 + j3 + 1) * f(j1 - m1) * f(j1 + m1) *
              f(j2 - m2) * f(j2 + m2)))

    S = 0.0
    for v in vmin:vmax
        S += (-1.0)^(v + j2 + m2) *
             (f(j2 + j3 + m1 - v) *
              f(j1 - m1 + v)) //
             (f(v) * f(j3 - j1 + j2 - v) * f(j3 + m3 - v) *
              f(v + j1 - j2 - m3))
    end
    return C * S
end

function clebsch_gordan(l1::Int, l2::Int, l3::Int)
    C = su2_clebsch_gordan(l1, l2, l3)
    Q1 = change_basis_real_to_complex(l1)
    Q2 = change_basis_real_to_complex(l2)
    Q3 = change_basis_real_to_complex(l3)
    result = zeros(ComplexF64, size(Q1, 2), size(Q2, 2), size(Q3, 2))
    Q3_conj_T = conj(transpose(Q3))
    @tensor result[j, l, m] := Q1[i, j] * Q2[k, l] * Q3_conj_T[m, n] * C[i, k, n]
    @assert all(abs.(imag.(C)) .< 1e-5)
    return real.(result)
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

Broadcast.broadcast(::typeof(wigner_D), l, args...) = broadcast(wigner_D, Ref(l), args...)
