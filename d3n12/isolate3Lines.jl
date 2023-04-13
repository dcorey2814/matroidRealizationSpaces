# ARGS[1] = input data file
# ARGS[2] = output data file

using Combinatorics
using Oscar

include("src/isolate3Lines.jl")

d3n12 = vec(readlines(joinpath(currentDir, ARGS[1])));

n3C12 = collect(powerset(1:12, 3, 3));
n3C12 = sort(n3C12, by =  x-> reverse(x));

io = open(ARGS[2], "w")


for c in cursor
    M = Matroid(c)
    ns = count_3_lines_thru_all_points(M)
    
    if length(ns) == 0
        continue
    end

    if minimum(ns) >= 3
        write(io, to_revlex(M, n3C12) * "\n")
    end
end

close(io)
