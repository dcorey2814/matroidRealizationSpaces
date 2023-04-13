# ARGS[1] = input data file
# ARGS[2] = output data file

using Combinatorics
using Oscar


currentDir = pwd()
include(joinpath(currentDir, "src/isolate3Lines.jl"));

d3n12 = vec(readlines(joinpath(currentDir, ARGS[1])));

n3C12 = collect(powerset(1:12, 3, 3));
n3C12 = sort(n3C12, by =  x-> reverse(x));

io = open(ARGS[2], "w")

i=0
for Mstr in d3n12
    M = matroid_from_revlex_basis_encoding(Mstr, 3, 12)    
    ns = count_3_lines_thru_all_points(M)
    
    if length(ns) == 0
        continue
    end

    if minimum(ns) >= 3
        write(io, to_revlex(M, n3C12) * "\n")
    end
    #global i += 1
    #println(i)
end

close(io)
