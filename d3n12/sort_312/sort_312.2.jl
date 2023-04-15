
using Combinatorics
using Oscar

#cd("..")

include("../src/matroid_realization.jl");
include("../src/reduceIdealsv2.jl");


d3n12_3_lines_2 = vec(readlines("matroids_d3n12/d3n12_3lines.2.dat"));

#checks for integer generator

function int_gens(I)
    ints = [x for x in I if is_constant(x)]
    return length(ints)>1    
end


#classification over QQ

#for z in 1:200
for z in 1:length(d3n12_3_lines_2)
    
    Mzstr = d3n12_3_lines_2[z]
    
    Mz = matroid_from_revlex_basis_encoding(Mzstr, 3, 12)
    I = matroid_to_reduced_expression(Mz, QQ)
    
    if I[1] == [0] #reduces to 0
        #println(z, " zero_ideal")
        
        open("matroids_d3n12/sorting2/zero_ideal.dat", "a")do file
        write(file, String(Mzstr),"\n")
        end
        continue
    
    elseif (1 in I[1] || -1 in I[1] || I == "Not Realizable 0 in Semigroup")
     #println(z, " not_realizable")   

        open("matroids_d3n12/sorting2/not_realizable.dat", "a")do file
        write(file, String(Mzstr),"\n")
        end
        continue
        
    elseif int_gens(I[1])
        #println(z, " int_generator")
    
        open("matroids_d3n12/sorting2/int_generator.dat", "a")do file
        write(file, String(Mzstr),"\n")
        end
        continue
        
    elseif length(I[1]) == 2
    	
        if (length(vars(I[1][2])) == 1 && degree(I[1][2],vars(I[1][2])[1]) >=2)
            #push!(simplified_principle_reducible,z)
            #println(z, " principal red")
            open("matroids_d3n12/sorting2/simplified_principle_reducible.dat", "a")do file
            write(file, String(Mzstr),"\n")
            end
            continue
                
        else
            #println(z, " principal n_red")
            #push!(simplified_principle_n_reducible,z)
            open("matroids_d3n12/sorting2/simplified_principle_n_reducible.dat", "a")do file
            write(file, String(Mzstr),"\n")
            end
            continue
        
        end
            
    elseif length(I[1])>2
    	#println(z, " not principal")
 
        open("matroids_d3n12/sorting2/simplified_not_principle.dat", "a")do file
        write(file, String(Mzstr),"\n")
        end
        continue
    else 
        #println(z, " net")
        #push!(net,z)
        open("matroids_d3n12/sorting2/net.dat", "a")do file
        write(file, String(Mzstr),"\n")
        end
        continue
    end
end
