using e3nn.o3

mutable struct IrrepsArray{T <: Real}
    irreps::Irreps
    array::Array{T}
    zero_flags::Union{Nothing, Tuple{Vararg{Bool}}}
    chunks::Union{Nothing, Vector{Union{Nothing, Array{T}}}}

    function IrrepsArray(irreps::Irreps, array::Array{T},
            zero_flags::Union{Nothing, Tuple{Vararg{Bool}}} = nothing,
            chunks::Union{Nothing, Vector} = nothing) where {T <: Real}
        # chunks=nothing) where {T<:Real}
        if size(array, ndims(array)) != dim(irreps)
            throw(ArgumentError("IrrepsArray: Array shape $(size(array)) incompatible with irreps $irreps"))
        end
        if !isnothing(chunks)
            for (mul_ir, chunk) in zip(irreps, chunks)
                if !isnothing(chunk) &&
                   size(chunk)[(end - 1):end] != (mul_ir.mul, dim(mul_ir.irrep))
                    throw(ArgumentError("IrrepsArray: chunk shape $(size(chunk)) incompatible with mul=$(mul_ir.mul) and dim(irrep)=$(dim(mul_ir.irrep))"))
                end
            end
        end
        if !isnothing(zero_flags) && length(zero_flags) != length(irreps)
            throw(ArgumentError("IrrepsArray: len(zero_flags) != len(irreps)"))
        end
        if !isnothing(chunks) && length(chunks) != length(irreps)
            throw(ArgumentError("IrrepsArray: len(chunks) != len(irreps)"))
        end
        new{T}(irreps, array,
            isnothing(zero_flags) ? repeat([false], length(irreps)) |> Tuple :
            Tuple(zero_flags),
            chunks)
    end
end

function IrrepsArray(irreps::String, array, zero_flags = nothing, chunks = nothing)
    IrrepsArray(Irreps(irreps), array, zero_flags, chunks)
end

function IrrepsArray(irreps, chunks, leading_shape)
    irreps = Irreps(irreps)
    if length(irreps) != length(chunks)
        throw(ArgumentError("IrrepsArray: len(irreps) != len(chunks)"))
    end
    # if any((x, mul_ir) -> !isnothing(x) && size(x) != (leading_shape..., (mul_ir.mul, dim(mul_ir.irrep))), zip(chunks, irreps))
    #     throw(ArgumentError("IrrepsArray: chunks shapes incompatible with leading shape and irreps"))
    # end
    if dim(irreps) > 0
        array = [if isnothing(x)
                     zeros(leading_shape..., dim(mul_ir))
                 else
                     reshape(x, leading_shape..., dim(mul_ir))
                 end
                 for (x, mul_ir) in zip(chunks, irreps)]
        array = cat(array..., dims = 1)
    else
        array = zeros(leading_shape..., 0)
    end
    zero_flags = (isnothing(x) for x in chunks) |> Tuple
    return IrrepsArray(irreps, array, zero_flags, chunks)
end

function IrrepsArray(arr::Array)
    if ndims(arr) == 0
        throw(ArgumentError("IrrepsArray: Cannot convert an array of rank 0 to an IrrepsArray."))
    end
    return IrrepsArray("$(size(arr, ndims(arr)))x0e", arr)
end

function chunks(ir_arr::IrrepsArray)
    if isnothing(ir_arr.chunks)
        leading_shape = size(ir_arr.array)[1:(end - 1)]
        zeros = ir_arr.zero_flags

        if length(ir_arr.irreps) == 1
            mul, ir = ir_arr.irreps[1]
            if zeros[1]
                return [nothing]
            else
                return [reshape(ir_arr.array, leading_shape..., mul, dim(ir))]
            end
        else
            # return [
            #     if z
            #         nothing
            #     else
            #         reshape(ir_arr.array[..., i], leading_shape..., mul, dim(ir))
            #     end
            #     for (z, i), (mul, ir) in zip(enumerate(ir_arr.irreps), ir_arr.irreps)
            # ]
            return []
        end
    end
end

Base.ndims(ir_arr::IrrepsArray) = ndims(ir_arr.array)
Base.length(ir_arr::IrrepsArray) = length(ir_arr.array)
Base.size(ir_arr::IrrepsArray) = size(ir_arr.array)

function Base.:+(x1::IrrepsArray, x2::Array)
    if all(mul_ir.irrep == "0e" for mul_ir in x1.irreps)
        return IrrepsArray(x1.irreps, x1.array .+ x2)
    else
        throw(ArgumentError("IrrepsArray($(x1.irreps), size=$(size(x1))) + scalar is not equivariant."))
    end
end

Base.:+(x1::IrrepsArray, x2::T) where {T <: Real} = x1 + [x2]

function Base.:+(x1::IrrepsArray, x2::IrrepsArray)
    if x1.irreps != x2.irreps
        throw(ArgumentError("IrrepsArray: Cannot add two IrrepsArrays with different irreps."))
    end
    zero_flags = [x && y for (x, y) in zip(x1.zero_flags, x2.zero_flags)] |> Tuple
    chunks = nothing
    if !isnothing(x1.chunks) && !isnothing(x2.chunks)
        chunks = [if isnothing(c1)
                      c2
                  else
                      if isnothing(c2)
                          c1
                      else
                          (c1 + c2)
                      end
                  end
                  for (c1, c2) in zip(x1.chunks, x2.chunks)]
    end
    return IrrepsArray(x1.irreps, x1.array + x2.array, zero_flags, chunks)
end

# TODO: needs chaning
function filter(keep, drop, lmax, ir_arr::IrrepsArray)
    if isnothing(keep) && isnothing(drop) && isnothing(lmax)
        return ir_arr
    end
    new_irreps = filter(ir_arr.irreps, keep, drop, lmax)
    return IrrepsArray(new_irreps, ir_arr.array, ir_arr.zero_flags, ir_arr.chunks)
end

const IntoIrreps = Union{Irrep, MulIrrep, String, Irreps}

function Base.show(io::IO, x::IrrepsArray)
    s = "$(x.irreps) $(x.array)"
    print(io, s)
end
