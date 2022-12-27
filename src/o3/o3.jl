module o3

include("irreps.jl")
export Irrep, Irreps, dim, isscalar, spherical_harmonics, simplify
export remove_zero_multiplicities, num_irreps, ls, lmax
# not implemented yet
export D_from_angles, D_from_quaternion, D_from_axis_angle, D_from_matrix

end