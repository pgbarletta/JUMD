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
    val = IntVector(undef, 0)

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
    val = IntVector(undef, 0)

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
    ids = Array{IntVector, 1}(undef, count)
    bot = 1
    for i = 1:count
        top = bot + val[i] - 1
        global ids[i] = indices[bot:top]
        bot += val[i]
    end

    return key, val, ids
end
"""
`weightedHist(v, bins, weights, density_bool, include_bounds_bool)`

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

julia> JUMD.weightedHist(v, bins, w)
  ([1.5, 2.5], [0.4, 0.6])
```
"""
function weightedHist(in_vec::RealVector, in_bins::RealVector,
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
`hisInd2D(x, y, bx, by, include_bounds_bool)`

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
function hisInd2D(in_vec_x::RealVector, in_vec_y::RealVector,
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

"""
`hisInd3D(x, y, z, bx, by, bz, include_bounds_bool)`

3D histogram of `x`, `y` and `z` given `bx`, `by` and `bz`  bins, but instead
of returning counts for each bin, it returns the indices of the elements that
matched for that particular bin, as a matrix of arrays.

If `include_bounds_bool == true` all ellements of `x`/`y`/`z` that fall beyond
the the boundary bins, will be assigned to to the corresponding boundary bin.

### Examples
```
TODO
```
"""
function hisInd3D(in_vec_x::RealVector, in_vec_y::RealVector, in_vec_z::RealVector,
    in_bins_x::RealVector, in_bins_y::RealVector, in_bins_z::RealVector,
    include_bounds::Bool = false)
    
    cnt = length(in_vec_x)
    if  length(in_vec_y) != cnt != length(in_vec_z)
        error("Input vectors length don't match. X: ", in_vec_x, " Y: ", in_vec_y,
        " Z: ", in_vec_z)
    end
         
    n_x = length(in_bins_x)
    n_y = length(in_bins_y)
    n_z = length(in_bins_z)
    
    his_ind = [Int64[] for i = 1:n_x, j = 1:n_y, z = 1:n_z]
    his = zeros(Int64, n_x, n_y, n_z)

    if include_bounds
        for i in 1:cnt
            x = searchsortedfirst(in_bins_x, in_vec_x[i])
            y = searchsortedfirst(in_bins_y, in_vec_y[i])
            z = searchsortedfirst(in_bins_z, in_vec_z[i])
            if x > n_x
                x = n_x
            end
            if y > n_y
                y = n_y
            end
            if z > n_z
                z = n_z
            end
            push!(his_ind[x, y, z], i)
            his[x, y, z] += 1
        end
    else
        for i in 1:cnt
            x = searchsortedfirst(in_bins_x, in_vec_x[i])
            y = searchsortedfirst(in_bins_y, in_vec_y[i])
            z = searchsortedfirst(in_bins_z, in_vec_z[i])
            if x > n_x || y > n_y || z > n_z
                continue
            end
            if (x == 1 && isless(x, in_vec_x[x])) || 
                (y == 1 && isless(y, in_vec_y[y])) ||
                (z == 1 && isless(z, in_vec_z[z]))
                continue
            end

            push!(his_ind[x, y, z], i)
            his[x, y, z] += 1
        end
    end
    
    return his_ind, his
end

"""
`suave(v, ws)`

Will smooth input vector `v` along a window of size `ws`. Resulting vector
will have type RealVector and length `floor(length(v) / ws)`. 

### Examples
```
julia> v = collect(1.:1.:10)
10-element Array{Float64,1}:
 1.0
 2.0
 3.0
 4.0
 5.0
 6.0
 7.0
 8.0
 9.0
10.0

julia> suave(v, 2)
5-element Array{Float64,1}:
 1.5
 3.5
 5.5
 7.5
 9.5
```
"""
function suave(v::RealVector, ws::Integer)::FloatVector
    vs = length(v) 
    ws_ = ws - 1
    return [ mean(v[i:i+ws_]) for i = 1:ws:length(v)-ws_ ]
end

"""
distancia(trj_fn::String, a::Integer, b::Integer, nchunks::Integer = 0, nframes::Integer = 0)

Read the trajectory at `trj_fn` using Chemfiles and returns a FloatVector
with the distance between atoms `a` and `b`. If `nchunks` != 0. the trajectory
will be divided in `nchunks` pieces to save RAM space. If `nframes` != 0, frames
1 to `nframes` will be read and the function will run slightly faster.

### Examples
```
NONE
```
"""
function distancia(trj_fn::String, a::Integer, b::Integer, nchunks::Integer = 0,
    nframes::Integer = 0)::FloatVector

    if nframes == 0
        in_trj = Trajectory(trj_fn)
        nframes = convert(Int64, nsteps(in_trj))
        close(in_trj)
    end
    distancias = Array{Float64, 1}(undef, nframes)

    if nchunks != 1
        st = convert(Int64, ceil(nframes / nchunks))
        # If, by any chance, `nframes` is divisble by `st` the last 2
        # elements will be nframes and the last iteration will be skipped.
        chunks = [ collect(0:st:nframes) ; nframes ]
    else
        chunks = [0; nframes]
    end

    for j in 1:nchunks
        in_trj = Trajectory(trj_fn)
        for i in chunks[j]:chunks[j+1]-1
            in_frm = read_step(in_trj, i)
            distancias[i+1] = norm(positions(in_frm)[:, a] - positions(in_frm)[:, b])
        end
        close(in_trj)
        GC.gc()
    end
    
    return distancias
end

"""
`triple(v::Array{Int64, 1})`

Will take vector of indices and triple them.

### Examples
```
julia> a = collect(3:5)
3-element Array{Int64,1}:
 3
 4
 5

julia> JUMD.triple(a)
9-element Array{Int64,1}:
  7
  8
  9
 10
 11
 12
 13
 14
 15
```
"""
function triple(v::IntVector)::IntVector
    n = length(v) 
    v_3 = IntVector(undef, n * 3)
    k = 1
    for i = 1:n
        v_3[k] = v[i] * 3 - 2
        v_3[k+1] = v[i] * 3 - 1 
        v_3[k+2] = v[i] * 3
        k+=3
    end
    return v_3
end
