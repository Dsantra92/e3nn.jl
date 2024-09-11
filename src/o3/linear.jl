struct Instruction
    i_in::Int
    i_out::Int
    path_shape::Tuple
    path_weight::Float64
    weight_std::Float64
end

struct Linear
    irreps_in::Irreps
    irreps_out::Irreps
    instructions::Vector{Instruction}
    output_mask::Vector{Bool}
end

# Default configuration dictionary
const __default_conf = Dict(
    "irrep_normalization" => "component",  # "component" or "norm"
    "path_normalization" => "element",  # "element" or "path"
    "gradient_normalization" => "path",  # "element", "path" or float,
    # "element" is the default in modules provided by pytorch/haiku
    "spherical_harmonics_algorithm" => "automatic",
    "spherical_harmonics_normalization" => "component",
    "custom_einsum_jvp" => false,
    "fused" => false,
    "sparse_tp" => false
)

__conf = deepcopy(__default_conf)

function config(name::String, value::Union{Nothing, Any} = nothing)
    if value === nothing
        return __conf[name]
    elseif haskey(__conf, name)
        __conf[name] = value
    else
        throw(ArgumentError("Unknown configuration option: $name"))
    end
end

function Linear(
        irreps_in::Irreps,
        irreps_out::Irreps;
        instructions::Union{Nothing, Vector{Tuple{Int, Int}}} = nothing,
        biases::Union{Nothing, Vector{Bool}, Bool} = nothing,
        path_normalization::Union{String, Float64, Nothing} = nothing,
        gradient_normalization::Union{String, Float64, Nothing} = nothing
)
    if isnothing(path_normalization)
        path_normalization = config("path_normalization")
    end
    if isa(path_normalization, String)
        path_normalization = Dict("element" => 0.0, "path" => 1.0)[path_normalization]
    end

    if isnothing(gradient_normalization)
        gradient_normalization = config("gradient_normalization")
    end
    if isa(gradient_normalization, String)
        gradient_normalization = Dict("element" => 0.0, "path" => 1.0)[gradient_normalization]
    end

    if isnothing(instructions)
        instructions = [(i_in, i_out) for (i_in, mul_ir_in) in enumerate(irreps_in)
                        for
                        (i_out, mul_ir_out) in enumerate(irreps_out)
                        if mul_ir_in == mul_ir_out]
    end

    instructions = [Instruction(
                        i_in,
                        i_out,
                        (irreps_in[i_in].mul, irreps_out[i_out].mul),
                        1.0, # dummy
                        1.0 # dummy
                    ) for (i_in, i_out) in instructions]

    function alpha(this)
        x = irreps_in[this.i_in].mul^path_normalization * sum(
            irreps_in[other.i_in].mul^(1.0 - path_normalization)
        for
        other in instructions if other.i_out == this.i_out
        )
        return x > 0 ? 1 / x : 1.0
    end

    instructions = [Instruction(
                        ins.i_in,
                        ins.i_out,
                        ins.path_shape,
                        sqrt(alpha(ins))^gradient_normalization,
                        sqrt(alpha(ins))^(1.0 - gradient_normalization)
                    ) for ins in instructions]
    if isnothing(biases)
        biases = fill(false, length(irreps_out))
    end
    if isa(biases, Bool)
        biases = [biases && isscalar(Irrep(ir)) for (_, ir) in irreps_out]
    end

    @assert length(biases) == length(irreps_out)
    @assert all(isscalar(mul_ir.irrep) || !b for (b, mul_ir) in zip(biases, irreps_out))

    append!(
        instructions,
        [Instruction(
             -1,
             i_out,
             (dim(mul_ir),),
             1.0,
             0.0
         ) for (i_out, (bias, mul_ir)) in enumerate(zip(biases, irreps_out)) if bias]
    )

    output_mask = if dim(irreps_out) > 0
        vcat(
            [if any(
                 (ins.i_out == i_out) && (0 ∉ ins.path_shape) for ins in instructions
             )
                 fill(true, dim(mul_ir))
             else
                 fill(false, dim(mul_ir))
             end
             for (i_out, mul_ir) in enumerate(irreps_out)]...,
        )
    else
        Bool[]
    end

    Linear(
        irreps_in,
        irreps_out,
        instructions,
        output_mask
    )
end

function Base.getproperty(fl::Linear, sym::Symbol)
    if sym === :num_weights
        return sum(prod(i.path_shape) for i in fl.instructions)
    else
        return getfield(fl, sym)
    end
end

# function (l::Linear)(ff

function Base.show(io::IO, fl::Linear)
    print(io,
        "Linear($(fl.irreps_in) -> $(fl.irreps_out), $(length(fl.instructions)) instructions, $(fl.num_weights) weights)")
end
