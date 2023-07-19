#resort 3-10 using using matroid_to_reduced_expression_min_chart   
#keep track of which chart we use

using Oscar
using Combinatorics
pm = Polymake


include("../../src/matroid_realization.jl");
include("../../src/reduceIdealsv2.jl")
include("../../src/JacobianCriterion.jl")


d3n10 = vec(readlines("../3lines_d3n10.dat"));

println("check function1")
function int_gens(I)
    ints = [x for x in I if is_constant(x)]
    return length(ints)>0    
end

println("check function 2")
function set_to_stars(set,n)
    
    z = []
    
    for a in 1:n
        if a in set
            push!(z,"*")
        else
            push!(z,"0")
        end
    end
    return join(z)
end


# open("2_to_recheck_nprinc_mult_singular_4_9.dat", "a") do file
               # write(file,data,"\n")
            #end

for z in 1:length(d3n10)
    println(z)
    Mz = matroid_from_revlex_basis_encoding(d3n10[z], 3, 10)
    I = matroid_to_reduced_expression_min_chart(Mz, QQ)
    data = join([set_to_stars(I[3],10),d3n10[z]])

 #check realizability   
    if I isa String    
        println("not realizable")
        open("not_realizable_3_10.dat", "a") do file
            write(file,data,"\n")
        end

    elseif isone(ideal(I[1])) 
               println("not realizable unit ideal")

               open("not_realizable_3_10.dat", "a") do file
                     write(file,data,"\n")
                end

    elseif (1 in I[1] || -1 in I[1])
            println("not realizable 1 or -1")
            open("not_realizable_3_10.dat", "a") do file
                write(file,data,"\n")
            end
        
#classify if realizable
    elseif (length(I[1]) == 0||I[1] == [0]) #reduces to 0 
            println("zero ideal")
            open("reduce_to_zero_3_10.dat", "a") do file
            write(file,data,"\n")
        end 
          

    elseif int_gens(I[1])#has integral generator. This only matters if reducing over ZZ
         println("int gen")
         open("int_generator_3_10.dat", "a") do file
            write(file,data,"\n")
        end
        
    elseif length(I[1]) == 1
        if (length(vars(I[1][1])) == 1)
             println("principal univatiate")#pricipal ideal in one variable
             open("principal_univariate_3_10.dat", "a") do file
                write(file,data,"\n")
             end


        else
            println("principal multivariate")#principle ideal in more than one variable
            open("principle_multivariate_3_10.dat", "a") do file
                write(file,data,"\n")
             end
     end   
    elseif length(I[1])>1
        
        Ivars = ideal_vars(I[1]) 
        if length(Ivars) == 1
                 println("not principal univariate")#univariate ideal with multiple generators
                 open("not_principle_univariate_3_10.dat", "a") do file
                     write(file,data,"\n")
                 end
        else
            println("not principle multivariate")#multivariate ideal with multiple generators
            open("not_principle_multivariate_3_10.dat", "a") do file
                    write(file,data,"\n")
         end
    end
    else
       println("to net")#to make sure sorting is exhaustive
       open("net_3_10.dat", "a") do file
            write(file,data,"\n")
     end
end
end

