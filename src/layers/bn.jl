using LinearAlgebra
using e3nn
using Statistics: mean, var

mutable struct BatchNorm{F, V, N, W}
    λ::F  # activation function
    β::V  # bias
    γ::V  # scale
    μ::W  # moving mean
    σ²::W # moving var
    ϵ::N
    momentum::N
    affine::Bool
    track_stats::Bool
    active::Union{Bool, Nothing}
    chs::Int # number of channels
    reduce::String
    normalization::String
    irreps::Irreps
end

# function batch_norm(
#     input::Irreps,
#     μ::Union{Nothing,AbstractVector},
#     σ²::Union{Nothing,AbstractVector},
#     ω::Union{Nothing,AbstractVector},
#     β::Union{Nothing,AbstractVector},
#     normalization::String,
#     reduce::String,
#     is_instance::Bool,
#     use_running_average::Bool,
#     use_affine::Bool,
#     momentum::Float64,
#     epsilon::Float64,
#     mask::Union{Nothing,AbstractVector}=nothing
# )
#     function _roll_avg(curr, update)
#         return (1 - momentum) * curr + momentum * update
#     end

#     batch, size... = size(input)[1:end-1]
#     input = reshape(input, (batch, prod(size), :))

#     if !is_instance
#         new_ra_means = []
#         new_ra_vars = []
#     end

#     new_chunks = []

#     i_wei = 1  # index for running_var and weight
#     i_rmu = 1  # index for running_mean
#     i_bia = 1  # index for bias

#     for (mul, ir) in zip(input.irreps.mul, input.irreps.ir)
#         chunk = input.chunks[ir]
#         if chunk === nothing
#             if !is_instance
#                 push!(new_ra_means, zeros(Float64, mul))
#             end
#             i_rmu += mul

#             if !is_instance
#                 push!(new_ra_vars, ones(Float64, mul))
#             end

#             if use_affine && ir.is_scalar
#                 i_bia += mul
#             end

#             push!(new_chunks, chunk)
#         else
#             if ir.is_scalar
#                 if is_instance
#                     mean = reshape(mean(chunk, dims=2), (batch, mul))
#                 else
#                     if mask === nothing
#                         mean = reshape(mean(chunk, dims=(1, 2)), mul)
#                     else
#                         mean = sum(mean(chunk, dims=2)[:, :] .* mask[:, :], dims=1) ./ sum(mask)
#                     end
#                     push!(new_ra_means, _roll_avg(ra_mean[i_rmu:i_rmu+mul-1], mean))
#                 end

#                 if use_running_average
#                     mean = ra_mean[i_rmu:i_rmu+mul-1]
#                 end
#                 i_rmu += mul

#                 chunk = chunk .- reshape(mean, (1, 1, mul, 1))
#             end

#             if normalization == "norm"
#                 norm_squared = sum(chunk .^ 2, dims=4)
#             elseif normalization == "component"
#                 norm_squared = mean(chunk .^ 2, dims=4)
#             else
#                 error("Invalid normalization option $normalization")
#             end

#             if reduce == "mean"
#                 norm_squared = mean(norm_squared, dims=2)
#             elseif reduce == "max"
#                 norm_squared = maximum(norm_squared, dims=2)
#             else
#                 error("Invalid reduce option $reduce")
#             end

#             if !is_instance
#                 if mask === nothing
#                     norm_squared = mean(norm_squared, dims=1)
#                 else
#                     norm_squared = sum(norm_squared .* mask[:, :], dims=1) ./ sum(mask)
#                 end
#                 push!(new_ra_vars, _roll_avg(ra_var[i_wei:i_wei+mul-1], norm_squared))
#             end

#             if use_running_average
#                 norm_squared = ra_var[i_wei:i_wei+mul-1]
#             end

#             inverse = 1 ./ sqrt.((1 - epsilon) .* norm_squared .+ epsilon)

#             if use_affine
#                 sub_weight = weight[i_wei:i_wei+mul-1]
#                 inverse = inverse .* sub_weight
#             end

#             chunk = chunk .* reshape(inverse, (1, 1, mul, 1))

#             if use_affine && ir.is_scalar
#                 sub_bias = bias[i_bia:i_bia+mul-1]
#                 chunk = chunk .+ reshape(sub_bias, (1, 1, mul, 1))
#                 i_bia += mul
#             end

#             push!(new_chunks, chunk)
#         end
#         i_wei += mul
#     end

#     @assert weight === nothing || i_wei == length(weight) + 1
#     @assert ra_var === nothing || i_wei == length(ra_var) + 1
#     @assert bias === nothing || i_bia == length(bias) + 1
#     @assert ra_mean === nothing || i_rmu == length(ra_mean) + 1

#     output = from_chunks(input.irreps, new_chunks, (batch, prod(size)), eltype(input))
#     output = reshape(output, (batch, size..., :))

#     if !is_instance
#         new_ra_means = !isempty(new_ra_means) ? vcat(new_ra_means...) : zeros(eltype(output), 0)
#         new_ra_vars = vcat(new_ra_vars...)
#     else
#         new_ra_means = nothing
#         new_ra_vars = nothing
#     end

#     return output, new_ra_means, new_ra_vars
# end
