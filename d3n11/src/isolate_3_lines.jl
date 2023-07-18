# Run from the top directory of this project

using Oscar 

currentDir = pwd() 
include(joinpath(currentDir, "src/fileHandling.jl"))
include(joinpath(currentDir, "src/isolate3Lines.jl"))
include(joinpath(currentDir, "src/matroid_realization.jl"))
include(joinpath(currentDir, "src/Jacobian_Criterion.jl"))


db = Polymake.Polydb.get_db();
collection = db["Matroids.Small"];
d3n11 = Polymake.Polydb.find(collection, Dict("RANK" => 3, "SIMPLE"=>true, "N_ELEMENTS"=>11));


n3C11 = subsets(collect(1:11), 3);
n3C11 = sort(n3C11, by =  x-> reverse(x));

lines_3 = []
for c in d3n11
    Q = Matroid(c)
    ns = count_3_lines_thru_all_points(Q)
    if length(ns) == 0
        continue
    end    
    if minimum(ns) >= 3
        push!(lines_3, to_revlex(Q, n3C11))
    end
end

lines_3_precomputed = vec(readlines(joinpath(currentDir, "d3n11/3lines_3_11.dat")));

print(Set(lines_3) == Set(lines_3_precomputed))
