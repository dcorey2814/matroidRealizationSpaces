using Oscar

currentDir = pwd() # make sure you are running from the main directory.
include(joinpath(currentDir, "src/fileHandling.jl"))
include(joinpath(currentDir, "src/matroid_realization.jl"))
include(joinpath(currentDir, "src/Jacobian_Criterion.jl"))

#pull database of simple 4-9 matroids
db = Polymake.Polydb.get_db();
collection = db["Matroids.Small"];
d4n9 = Polymake.Polydb.find(collection, Dict("RANK" => 4, "SIMPLE"=>true, "N_ELEMENTS"=>9));

#counts number of planes an element of the ground set is in
function in_planes(n,L)

    Ln = [l for l in L if n in l]

    return(length(Ln))

end



d4C9 = subsets(collect(1:9), 4);
d4C9 = sort(d4C9, by =  x-> reverse(x));

#isolate simple connected matroids satisfying 4 lines property
planes_4 = []

for t in d4n9
    
   # println(t)

    Qt = Matroid(t)

    L = [h for h in hyperplanes(Qt) if length(h)>3]

    ns = [n for n in 1:9 if in_planes(n,L)>3]

    if (length(ns) == 9 && is_connected(Qt))
                push!(planes_4,to_revlex(Qt,d4C9))
    end
end


four_planes_precomputed = vec(readlines("d4n9/4planes_connected_4_9.dat"))

print(Set(planes_4) == Set(four_planes_precomputed))
