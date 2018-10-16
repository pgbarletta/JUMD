module JUMD
using FileIO, Chemfiles, DelimitedFiles, LinearAlgebra
"""

JUMD.contar(v)


Devuelve los valores únicos dentro del **v**, junto a un conteo de estos.

### Examples
```
julia> veca = [ repeat(collect(0:3), 3); collect(4:5); repeat(collect(6:9), 5) ];
julia> JUMD.contar(veca)
([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [3, 3, 3, 3, 1, 1, 5, 5, 5, 5])

julia> Dict(zip(JUMD.contar(veca)...))
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

JUMD.contarIndexar(v)


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

# Agarra una matriz de PCA en Calpha de 3Nx3N-6 y devuelve una lista de 3N-6
# matrices, c/u de 3x3N. C/ mtx es 1 modo reordenado p/ matchear las matrices
# de coordenadas de carbonos alfa.
function format_pca_aa(in_mtx::Array{Float64, 2})
    m, n = size(in_mtx)
    aa = Int64
    try
        aa = convert(Int64, m / 3)
    catch
        error("Vector length: ", m, " is not divisible by 3.")
    end

    list_out_mtx = Array{Array{Float64, 2}, 1}(undef, n);
    [ list_out_mtx[j] = reshape(in_mtx[:, j], 3, aa) for j = 1:n ]
    
    return list_out_mtx
end

# Agarra la topología una matriz de PCA en Calpha de 3Nx3N-6 y devuelve una
# lista de 3N-6 matrices, c/u de 3xNatomos. C/ mtx es 1 modo reordenado p/
# matchear las matrices de coordenadas del pdb q dió lugar a la topología.
# También devuelve un array con el nro de atomos q tiene c/ aa
function format_pca_atom(in_top::Topology, in_mtx::Array{Float64, 2},
    mask::Array{Float64, 1} = 0)
    # Preparo variables
    aa = Int64
    aa_3 = Int64
    if mask == 0 
        aa = convert(Int64, count_residues(in_top))
    else
        aa = length(mask)
    end
    aa_3 = aa * 3
    
    m, n = size(in_mtx)
    if m != aa_3
        error("Input vector with wrong dimensions: ", m, "  ", (aa_3, 1))
    end

    # Determino orden de residuos (hay q actualizar el Julia Chemfiles)
    tmp = Array{Int64}(undef, aa)
    ids = Array{Int64}(undef, aa)
    [ ids[i+1] = convert(Int64, id((Residue(in_top, i)))) for i = 0:aa-1 ]
    idx = sortperm(ids)
    # Determino el nro de atomos de c/ aminoácido. Resto 1 pq Chemfiles tiene 0-indexing
    [ tmp[i] = size(Residue(in_top, mask[i] - 1)) for i = 1:aa ]
    natom_aa = tmp[idx]
    natoms = sum(natom_aa)

    # Adapto el vector p/ darle la misma forma q la matriz de coordenadas
    list_out_mtx = Array{Array{Float64, 2}, 1}(n);
    
    for j in 1:n
        vector = reshape(in_mtx[:, j], 3, aa)
        list_out_mtx[j] = Array{Float64}(3, natoms)
        cursor = 0
        for i = 1:aa
            rango = Array{Int64}(natom_aa[i])
            if i == 1
                list_out_mtx[j][:, 1:natom_aa[i]] = repmat(vector[:, i], 1, natom_aa[i])
                cursor = natom_aa[i]
                continue
            end
            rango = collect(cursor+1:cursor + natom_aa[i])
            list_out_mtx[j][:, rango] = repmat(vector[:, i], 1, natom_aa[i])
            cursor += natom_aa[i]
        end
    end

    return list_out_mtx, natom_aa
end

function get_κ(in_vec::Array{Float64, 1})
    not_null = copy(in_vec)
    not_null[not_null .== 0] .= 0.000001
    κ = (exp.(-mapslices(x -> sum(x), mapslices(x->x.^2 .* log.(x.^2), not_null, dims = 1), dims = 1))
        / length(not_null))[1]
    return κ
end

function get_pnum(in_vec::Array{Float64, 1})
    nor_vec = in_vec ./ norm(in_vec) 
    return convert(Int64, round(sum(nor_vec .^ 4) .^ -1))
end

function tognm(vtor_anm::Array{Float64, 1})
    vtor_gnm = Array{Float64, 1}
    try
        vtor_gnm = Array{Float64, 1}(undef, convert(Int64, length(vtor_anm)/3));
    catch e
        warn("Input vector's length is not a 3 multiplier")
        error(e)
    end
    vtor_anm =  vtor_anm.^2
    n = convert(Int64, length(vtor_anm)/3)
    
    [ vtor_gnm[i] = sqrt(vtor_anm[i*3-2] + vtor_anm[i*3-1] + vtor_anm[i*3])
        for i = 1:n ]
    
    return vtor_gnm
end

function WeightedHist(in_vec, in_bins, in_weight, density = false, include_bounds = true)
    # Safety check    
    if length(in_vec) != length(in_weight)
        error("Each element of the input vector needs one weight")
        return
    end
    
    # Prepare variables
    nbins = length(in_bins) - 1
    out_counts = Array{Float64}(undef, nbins)
    
    # Get weighted histogram
    if include_bounds
        for i=1:nbins
            if i == 1
                # Include those that fall before the beggining of the bins
                temp_bool = (in_vec .>= in_bins[i]) .& (in_vec .< in_bins[i+1]) .| (in_vec .<= in_bins[i])
                out_counts[i] = sum(in_weight[temp_bool])
            elseif i == length(in_bins)-1
                # Include those that fall after the end of the bins
                temp_bool = (in_vec .>= in_bins[i]) .& (in_vec .< in_bins[i+1]) .| (in_vec .>= in_bins[end])
                out_counts[i] = sum(in_weight[temp_bool])
            else
                temp_bool = (in_vec .>= in_bins[i]) .& (in_vec .< in_bins[i+1])
                out_counts[i] = sum(in_weight[temp_bool])
            end
        end
    else
        for i=1:length(in_bins)-1
            temp_bool = (in_vec .>= in_bins[i]) .& (in_vec .< in_bins[i+1])
            out_counts[i] = sum(in_weight[temp_bool])
        end
    end
    
    # Get bins middle points
    out_middle = (in_bins[1:end-1] + in_bins[2:end]) / 2
    
    # Turn counts into density
    if (density == true)
        out_counts = out_counts ./ sum(out_counts) 
    end
    return out_middle, out_counts
end

function MatHisInd2D(in_vec_x::AbstractArray, in_vec_y::AbstractArray,
    in_bins_x::AbstractArray, in_bins_y::AbstractArray,
    include_bounds = true)

    cnt = length(in_vec_x)
    if  length(in_vec_y) != cnt
        error("Input vectors length don't match. X: ", in_vec_x, " Y: ", in_vec_y)
    end
         
    n_x = length(in_bins_x)
    n_y = length(in_bins_y)
    
    his_ind = [Int[] for i=1:n_x, j=1:n_y]
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

#
function read_ptraj_modes(filename, nmodes::Int64 = 0, norma::Bool = true)
    modes_text = readdlm(filename, skipstart=0, skipblanks=true, comments=true,
        comment_char='*')

    if nmodes == 0
        nmodes = convert(Int64, modes_text[1, 5])
    end
    ncoords = convert(Int64, modes_text[2, 1])

    lines = ceil(Int64, ncoords/7)
    fields_per_line = lines * 7
    evalue = Array{Float64, 1}(undef, nmodes)
    mode = Array{Float64, 2}(undef, ncoords, nmodes)
    
    j = lines + 1 + 2 # 1 p/ q lea la prox linea 2 por el header
    for i=1:nmodes
        evalue[i] = modes_text[j, 2]
        temp = permutedims(modes_text[(j+1):(lines+j), :], [2, 1])
        temp2 = reshape(temp, fields_per_line)
        [ pop!(temp2) for k = ncoords+1:fields_per_line ]
        mode[:, i] = convert(Array{Float64, 1}, temp2)
        j = j + lines + 1
    end

    if norma == true
        for i=1:nmodes
            mode[: ,i] = mode[:, i] / norm(mode[:, i])
        end
    end

     return mode, evalue
end

function energia_gdte(evals::Array{Float64, 1}, gdte::Array{Float64, 1}, d::Float64 = 1.)

    if length(evals) != length(gdte)
        error("Lengths of evals and gdte don't match. Aborting.") 
    end
    
    # Declaro cte de boltzmann, avogadro, y temperatura.
    k = 1.38064852e-23
    avgdro = 6.0221409e+23
    T = 298
    RT =  k * avgdro * T * 1E-3 * 0.239006 # Kcal/mol
    cte = 11792.08316093831
    
    return d^2 * 0.5 * RT * sum(evals.^2 .* gdte.^2) / cte # Kcal/mol
end

end
