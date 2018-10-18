"""
`contar(v)`

Devuelve los valores únicos dentro del `v`, junto a un conteo de estos.

### Examples
```
julia> veca = [ repeat(collect(0:3), 3); collect(4:5); repeat(collect(6:9), 5) ];
julia> contar(veca)
([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [3, 3, 3, 3, 1, 1, 5, 5, 5, 5])

julia> Dict(zip(contar(veca)...))
Dict{Int64,Int64} with 10 entries:
  0 => 3
  4 => 1
  7 => 5
  9 => 5
  2 => 3
  3 => 3
  5 => 1
  8 => 5
  6 => 5
  1 => 3
```
"""
function contar(in_vtor::Array{T, 1}) where T
    in_v = sort(in_vtor)
    n = length(in_v)
    key = Array{T, 1}(undef, 0)
    val = Array{Int64, 1}(undef, 0)

    let st = 1
        for k = 1:n
            temp = in_v[st:n]
            idx = searchsortedlast(temp, in_v[st])
            st += idx
        
            push!(key, temp[idx])
            push!(val, idx)
            st > n && break
        end
    end

    return key, val
end

"""
`JUMD.contarIndexar(v)`

Devuelve los valores únicos dentro del v, junto a un conteo de estos y
un array con los indices donde aparece c/ valor único

### Examples
```
julia> veca = [ repeat(collect(0:3), 3); collect(4:5); repeat(collect(6:9), 5) ];
julia> nombre, frecu, indices = contarIndexar(veca)
([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [3, 3, 3, 3, 1, 1, 5, 5, 5, 5],Array{Int64,1}[[1, 5, 9], [2, 6, 10], [3, 7, 11], [4, 8, 12], [13], [14], [15, 19, 23, 27, 31], [16, 20, 24, 28, 32], [17, 21, 25, 29, 33], [18, 22, 26, 30, 34]])

julia> veca[indices[1]]
3-element Array{Int64,1}:
 0
 0
 0
julia> veca[indices[5]]
1-element Array{Int64,1}:
 4
```
"""
function contarIndexar(in_vtor::Array{T, 1}) where T
    indices = sortperm(in_vtor)
    in_v = in_vtor[indices]
    n = length(in_v)
    key = Array{T, 1}(undef, 0)
    val = Array{Int64, 1}(undef, 0)

    let st = 1
        for k = 1:n
            temp = in_v[st:n]
            idx = searchsortedlast(temp, in_v[st])
            st += idx

            push!(key, temp[idx])
            push!(val, idx)
            st > n && break
        end
    end
    
    count = length(val)
    ids = Array{Array{Int64, 1}, 1}(undef, count)
    bot = 1
    for i = 1:count
        top = bot + val[i] - 1
        global ids[i] = indices[bot:top]
        bot += val[i]
    end

    return key, val, ids
end
"""
`WeightedHist(v, bins, weights, density_bool, include_bounds_bool)`

Histogram of `v`, given bins `bins` but instead of counting each element in
`v` as 1, they are counted as their corresponding weight in `weights`.
The sum of `weights` vector should amount to 1.

If `density_bool == true`, then the sum of all counts will amount to 1.

If `include_bounds_bool == true` all ellements of `v` smaller(larger) than
the first(last) bin, will be assigned to it.

### Examples
```
julia> v = [1.1; 1.2; 1.3; 2.1; 2.2]
5-element Array{Float64,1}:
 1.1
 1.2
 1.3
 2.1
 2.2

julia> w = [.1; .15 ; .15; .3 ; .3]
5-element Array{Float64,1}:
  0.1 
  0.15
  0.15
  0.3 
  0.3 
 
julia> bins = [1.; 2.; 3.]
3-element Array{Float64,1}:
   1.0
   2.0
   3.0

julia> JUMD.WeightedHist(v, bins, w)
  ([1.5, 2.5], [0.4, 0.6])
```
"""
function WeightedHist(in_vec::RealVector, in_bins::RealVector,
    in_weight::RealVector, density::Bool = true, include_bounds::Bool = false)
    
    if length(in_vec) != length(in_weight)
        error("Each element of the input vector needs one weight")
        return
    end
    
    nbins = length(in_bins) - 1
    out_counts = Array{Float64, 1}(undef, nbins)
    
    # Get weighted histogram
    if include_bounds
        for i=1:nbins
            if i == 1
                # Include those that fall before the beggining of the bins
                temp_bool = (in_vec .>= in_bins[i]) .& 
                    (in_vec .< in_bins[i+1]) .| (in_vec .<= in_bins[i])
                out_counts[i] = sum(in_weight[temp_bool])
            elseif i == nbins
                # Include those that fall after the end of the bins
                temp_bool = (in_vec .>= in_bins[i]) .&
                    (in_vec .< in_bins[i+1]) .| (in_vec .>= in_bins[end])
                out_counts[i] = sum(in_weight[temp_bool])
            else
                temp_bool = (in_vec .>= in_bins[i]) .&
                    (in_vec .< in_bins[i+1])
                out_counts[i] = sum(in_weight[temp_bool])
            end
        end
    else
        for i=1:nbins
            temp_bool = (in_vec .>= in_bins[i]) .& (in_vec .< in_bins[i+1])
            out_counts[i] = sum(in_weight[temp_bool])
        end
    end
    
    # Get bins middle points.
    out_middle = (in_bins[1:end-1] + in_bins[2:end]) ./ 2
    
    # Turn counts into density, if required.
    if (density)
        out_counts = out_counts ./ sum(out_counts) 
    end
    
    return out_middle, out_counts
end
"""
`HisInd2D(x, y, bx, by, include_bounds_bool)`

2D histogram of `x` and `y` given `bx` and `by` bins, but instead of returning
counts for each bin, it returns the indices of the elements that matched for
that particular bin, as a matrix of arrays.

If `include_bounds_bool == true` all ellements of `x`/`y` that fall beyond the
the boundary bins, will be assigned to to the corresponding boundary bin.

### Examples
```
TODO
```
"""
function HisInd2D(in_vec_x::RealVector, in_vec_y::RealVector,
    in_bins_x::RealVector, in_bins_y::RealVector, include_bounds::Bool = false)
    
    cnt = length(in_vec_x)
    if  length(in_vec_y) != cnt
        error("Input vectors length don't match. X: ", in_vec_x, " Y: ", in_vec_y)
    end
         
    n_x = length(in_bins_x)
    n_y = length(in_bins_y)
    
    his_ind = [Int64[] for i=1:n_x, j=1:n_y]
    his = zeros(Int64, n_x, n_y)
    if include_bounds
        for i in 1:cnt
            x = searchsortedfirst(in_bins_x, in_vec_x[i])
            y = searchsortedfirst(in_bins_y, in_vec_y[i])
            if x > n_x
                x = n_x
            end
            if y > n_y
                y = n_y
            end
            
            push!(his_ind[x, y], i)
            his[x, y] += 1
        end
    else
        for i in 1:cnt
            x = searchsortedfirst(in_bins_x, in_vec_x[i])
            y = searchsortedfirst(in_bins_y, in_vec_y[i])
        
            if x > n_x || y > n_y
                continue
            end
            if (x == 1 && isless(x, in_vec_x[x])) || (y == 1 && isless(y, in_vec_y[y]))
                continue
            end

            push!(his_ind[x, y], i)
            his[x, y] += 1
        end
    end
    
    return his_ind, his
end