using Combinatorics
using Oscar

currentDir = pwd()

include(joinpath(currentDir, "src/isolate3Lines.jl"))

db = Polymake.Polydb.get_db();
collection = db["Matroids.Small"];
cursor=Polymake.Polydb.find(collection, Dict("RANK" => 3, "SIMPLE"=>true, "N_ELEMENTS"=>11));


n3C11 = collect(powerset(1:11, 3, 3));
n3C11 = sort(n3C11, by =  x-> reverse(x));

io = open(joinpath(currentDir, "d3n11/3lines_d3n11.dat"), "w")

for c in cursor
    M = Matroid(c)
    ns = count_3_lines_thru_all_points(M)
    
    if length(ns) == 0
        continue
    end

    if minimum(ns) >= 3
        write(io, to_revlex(M, n3C11) * "\n")
    end

end

close(io)
