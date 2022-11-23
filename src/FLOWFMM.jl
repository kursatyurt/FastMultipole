module FLOWFMM

import Base.:^
# import Statistics as S
using LinearAlgebra
using StaticArrays

# const THETA = 4
const ONE_OVER_4PI = 1/4/pi

# indices
const i_POSITION = 1:3 # inside .bodies member
const i_SIGMA = 4 # inside .bodies member; body radius
const i_POTENTIAL = 1:4 # inside .potential member
const i_POTENTIAL_JACOBIAN = 5:16 # inside .potential member
const i_POTENTIAL_HESSIAN = 17:52 # inside .potential member

for file in ["misc", "options", "derivatives", "element", "tree", "direct", "spherical", "fmm"]
    include(file*".jl")
end

end
