using Oscar
using Combinatorics

function lines(M::Matroid)
    return [l for l in hyperplanes(M) if length(l) > rank(M)-1 ]
end

function count_lines_thru_one_point(Ls,i)
    return length([l for l in Ls if i in l])
end

function count_3_lines_thru_all_points(M)
    Ls = lines(M)    
    return [count_lines_thru_one_point(Ls,i) for i in matroid_groundset(M)]
end



function to_revlex(M, nCd)
    B = bases(M)
    l = []
    for b in nCd
        if b in B
            push!(l,"*")
        else
            push!(l,"0")
        end
    end
    return join(l)
end



db = Polymake.Polydb.get_db();
collection = db["Matroids.Small"];
cursor=Polymake.Polydb.find(collection, Dict("RANK" => 3,"SIMPLE"=>true,"N_ELEMENTS"=>11));

n3C11 = collect(powerset(1:11, 3, 3));
n3C11 = sort(n3C11, by =  x-> reverse(x));

io = open("d3n11_3lines.dat", "w")

i = 0;
for c in cursor
    M = Matroid(c)
    ns = count_3_lines_thru_all_points(M)
    
    if length(ns) == 0
        continue
    end

    if minimum(ns) >= 3
        write(io, to_revlex(M, n3C11) * "\n")
    end

    # checkpoint
    if (i%5000 == 0)
        close(io)
        io = open("d3n11_3lines.dat", "a")
    end

    i = i+1
#    println(i)
end


close(io)
