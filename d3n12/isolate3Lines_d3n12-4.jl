using Combinatorics
using Oscar

include("../src/isolate3Lines.jl")

db = Polymake.Polydb.get_db();
collection = db["Matroids.Small"];

#cursor=Polymake.Polydb.find(collection, Dict("RANK" => 3,"SIMPLE"=>true,"N_ELEMENTS"=>12));

cursor=Polymake.Polydb.find(collection, Dict("RANK" => 3,"SIMPLE"=>true,"N_ELEMENTS"=>12), opts=Dict("skip"=>12000001, "limit"=>16000000));

n3C12 = collect(powerset(1:12, 3, 3));
n3C12 = sort(n3C12, by =  x-> reverse(x));

io = open("matroids_d3n12/d3n12_3lines.4.dat", "w")

#i = 0;
for c in cursor
    M = Matroid(c)
    ns = count_3_lines_thru_all_points(M)
    
    if length(ns) == 0
        continue
    end

    if minimum(ns) >= 3
        write(io, to_revlex(M, n3C12) * "\n")
    end

    # checkpoint
#    if (i%100000 == 0)
#        close(io)
#        global io = open("matroids_d3n12/d3n12_3lines.4.dat", "a")
#        println(i)
#    end

#    global i += 1

end

close(io)
