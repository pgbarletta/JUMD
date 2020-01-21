using Chemfiles
using LinearAlgebra, Statistics
using FileIO, Printf

# nombres de archivos de salida y número de frames
evecs_fn = "evecs"
evals_fn = "evals"
K = 598

# Arranco con el primer frame
xyz = Array{Array{Float64, 2}, 1}(undef, K)
tr = Chemfiles.Trajectory("snapshot1.pdb")

fr = read(tr)
top = Topology(fr)
n = convert(Int64, size(top))
cas = Array{Int64, 1}()
for i in 0:n-1
    if name(Atom(top, i)) == "CA"
        push!(cas, i+1)
    end
end
aa = length(cas)
aa3 = 3*aa

tmp = positions(fr)[:, cas]
close(tr)
# xyz tiene 3 columnas p/ c/ CA (X, Y y Z). C/ fila son las coordenadas de los
# CAs en 1 snapshot.
xyz = Array{Float64, 2}(undef, K, aa3)
for i in eachindex(tmp)
     xyz[1, i] = tmp[i]
end

# Ahora agrego las coordenadas de c/ frame a xyz
for i in 2:K
    tr = Chemfiles.Trajectory(string("snapshot", i, ".pdb"))
    fr = read(tr)
    tmp = positions(fr)[:, cas]
    for k in eachindex(tmp)
        xyz[i, k] = tmp[k]
    end
    close(tr)
end

# Calculo matriz de covarianza y luego la diagonalizo. Finalmente, descarto
# desplazamiento y rotación
covar = [ cor(xyz[:, i], xyz[:, j]) for i in 1:aa3, j in 1:aa3  ]
evals, evecs = eigen(covar);
evals = evals[7:end]
evecs = evecs[:, 7:end];

# Escribo salida con buen formato
io = open(evecs_fn, "w")
for i in 1:aa3
    [ @printf(io, "%8.4f ", evecs[i, j]) for j in 1:aa3-6 ]
    @printf(io, "\n")
end
close(io)

io = open(evals_fn, "w")
[ @printf(io, "%14.9f\n", evals[j]) for j in 1:aa3-6 ]
close(io)
