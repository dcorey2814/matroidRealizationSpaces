#println(pwd())

using Combinatorics
using Oscar
pm = Polymake

#cd("..")

include("src/matroid_realization.jl");
include("src/reduceIdealsv2.jl");

d3n11 = vec(readlines("3_11/d3n11.dat"));
d3n11_3_lines = vec(readlines("3_11/d3n11_3lines.dat"));

#checks for integer generator

function int_gens(I)
    ints = [x for x in I if is_constant(x)]
    return length(ints)>1    
end


#classification over QQ
for z in 1:length(d3n11_3_lines)
    #println(z)
    
    #println(z,"\n\n")
  
    Mz = matroid_from_revlex_basis_encoding(d3n11_3_lines[z], 3, 11)
    I = matroid_to_reduced_expression(Mz, QQ)
    #println(z, " success")
    
    #println("ch1")
    if I[1] == [0] #reduces to 0
        println(z, " zero_ideal")
        #push!(zero_ideal,z)
        
        open("3_11/sorting/zero_ideal2.dat", "a")do file
        write(file, String(d3n11_3_lines[z]),"\n")
        end
        continue
    
    elseif (1 in I[1] || -1 in I[1] || I == "Not Realizable 0 in Semigroup")
     
     println(z, " not_realizable")   
    #println("ch2")
        
        #push!(not_realizable,z)
        open("3_11/sorting/not_realizable2.dat", "a")do file
        write(file, String(d3n11_3_lines[z]),"\n")
        end
        continue
        
    elseif int_gens(I[1])
        println(z, " int_generator")
               
    #println("ch3")
        
       # push!(int_generator,z)
    
        open("3_11/sorting/int_generator2.dat", "a")do file
        write(file, String(d3n11_3_lines[z]),"\n")
        end
        continue
        
    elseif length(I[1]) == 2
    	
            
        if (length(vars(I[1][2])) == 1 && degree(I[1][2],vars(I[1][2])[1]) >=2)
                
            #push!(simplified_principle_reducible,z)
            println(z, " principal red")
        open("3_11/sorting/simplified_principle_reducible2.dat", "a")do file
        write(file, String(d3n11_3_lines[z]),"\n")
        end
        continue
                
        else
                
        println(z, " principal n_red")
        #push!(simplified_principle_n_reducible,z)
        open("3_11/sorting/simplified_principle_n_reducible2.dat", "a")do file
        write(file, String(d3n11_3_lines[z]),"\n")
        end
        continue
        
        end
            
    elseif length(I[1])>2
    
    	println(z, " not principal")
        
            #push!(simplified_not_principle,z)
        
    	
        open("3_11/sorting/simplified_not_principle2.dat", "a")do file
        write(file, String(d3n11_3_lines[z]),"\n")
        end
        continue
    else
        
        println(z, " net")
        #push!(net,z)
        open("3_11/sorting/net2.dat", "a")do file
        write(file, String(d3n11_3_lines[z]),"\n")
        end
        continue
        
    end
end
