"""
`formatPcaAA(mtx)`

Takes a 3Nx3N-6 normal modes matrix (where N = nbr of particles) and gives a
3N-6 array of 3x3N matrices, each corresponding to a normal mode.

### Examples
```
TODO
```
"""
function formatPcaAA(in_mtx::RealMatrix)
    m, n = size(in_mtx)
    aa = Int64
    try
        aa = convert(Int64, m / 3)
    catch
        error("Vector length: ", m, " is not divisible by 3.")
    end

    list_out_mtx = Array{Array{T, 2}, 1}(undef, n);
    [ list_out_mtx[j] = reshape(in_mtx[:, j], 3, aa) for j = 1:n ]
    
    return list_out_mtx
end
"""
`formatPcaAtom(topology, mtx, mask)`

TODO
### Examples
```
TODO
```
"""
function formatPcaAtom(in_top::Topology, in_mtx::Array{Float64, 2},
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
"""
`getκ(v)`

Calculate the collectivity of the input normal mode `v` according to Sanjeouand.
[citation needed]

### Examples
```
TODO
```
"""
function getκ(in_vec::RealVector)
    not_null = copy(in_vec)
    not_null[not_null .== 0] .= 0.000001
    κ = (exp.(-mapslices(x -> sum(x), mapslices(x->x.^2 .* log.(x.^2), not_null, dims = 1), dims = 1))
        / length(not_null))[1]
    return κ
end
"""
`getPnum(v)`

Calculate the participation number of the input mode/difference vector
according to ... [citation needed]

### Examples

```
TODO
```

"""
function getPnum(in_vec::RealVector)
    nor_vec = in_vec ./ norm(in_vec) 
    return convert(Int64, round(sum(nor_vec .^ 4) .^ -1))
end
"""
`toGnm(v)`

Turn the cartesian (XYZ) mode into a Gaussian Normal Mode-like vector. That
is, turn xyz displacements into absolute displacements.

### Examples
```
TODO
```
"""
function toGnm(vtor_anm::RealVector)
    m = length(vtor_anm)
    n = Int64
    try
        n = convert(Int64, m/3)
    catch e
        error("Input vector's length is not divisible by 3.")
    end
    
    vtor_gnm = Array{Float64, 1}(undef, n);
        
    [ vtor_gnm[i] = sqrt(vtor_anm[i*3-2]^2 + vtor_anm[i*3-1]^2 + vtor_anm[i*3]^2)
        for i = 1:n ]
    
    return vtor_gnm
end
"""
`readPtrajModes(f, n, normalize_bool)`

Read the first `n` PCA modes from the `f` file. Assumes Amber format.
If `n == 0`, read all vectors.
If `normalize_bool == true` will normalize vector norm to 1.

### Examples
```
TODO
```
"""
function readPtrajModes(filename::AbstractString, nmodes::Int64 = 0,
    norma::Bool = true)
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

    if norma
        [ mode[: ,i] = mode[:, i] ./ norm(mode[:, i]) for i=1:nmodes ]
    end

     return mode, evalue
end
"""
`energiaGdte(e, g, d)`


Calculate the energy of the Volume Gradient Vector `g` given the modes
eigenvalues `e` for a displacement of size `d`.

### Examples
```
TODO
```
"""
function energiaGdte(evals::RealVector, gdte::RealVector, d::Float64 = 1.)

    if length(evals) != length(gdte)
        error("Lengths of evals and gdte don't match. Aborting.") 
    end
    
    # Declaro cte de boltzmann, avogadro, y temperatura.
    k = 1.38064852e-23
    avgdro = 6.0221409e+23
    T = 298
    RT =  k * avgdro * T * 1E-3 * 0.239006 # Kcal/mol
    cte = 11792.08316093831

    a = (RT/cte) * sum(evals.^2 .* gdte.^2)
    
    U = (1/2) * a * d^2 # Kcal/mol
    return U
end
"""
`energiaGdteNMA(e, g, d)`


Calculate the energy of the Volume Gradient Vector `g` given the modes
eigenvalues `e` for a displacement of size `d`. Eigenvalues should be
angular frequencies

### Examples
```
TODO
```
"""
function energiaGdteNMA(evals::RealVector, gdte::RealVector, d::Float64 = 1.)

    if length(evals) != length(gdte)
        error("Lengths of evals and gdte don't match. Aborting.") 
    end
    
    # Declaro cte de boltzmann, avogadro, y temperatura.
    
    U = d^2 * 0.5 * sum(evals.^2 .* gdte.^2) 
    return U
end
