module O3

include("irreps.jl")
export Irrep, Irreps
export dim, isscalar, lmax, ls, count, mul_gcd, num_irreps,
       remove_zero_multiplicities, simplify, sort, unify, regroup

end
