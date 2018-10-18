module JUMD
using FileIO, DelimitedFiles, LinearAlgebra, Chemfiles

const RealVector{T<:Real} =  Array{T, 1}
const RealMatrix{T<:Real} =  Array{T, 2}


include("utils.jl")
include("modes.jl")

end
