module JUMD
using FileIO, DelimitedFiles, LinearAlgebra, Statistics, Chemfiles

const RealVector{T<:Real} =  Array{T, 1}
const FloatVector{T<:Real} =  Array{Float64, 1}
const IntVector{T<:Real} =  Array{Int64, 1}
const RealMatrix{T<:Real} =  Array{T, 2}


include("utils.jl")
include("modes.jl")

end
