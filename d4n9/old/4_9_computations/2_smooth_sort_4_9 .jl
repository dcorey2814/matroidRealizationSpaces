#resort 4-9 using using matroid_to_reduced_expression_min_chart   
#keep track of which chart we use

using Oscar
using Combinatorics
pm = Polymake


include("../../src/matroid_realization.jl");
include("../../src/reduceIdealsv2.jl")
include("../../src/JacobianCriterion.jl")


d4n9 = vec(readlines("../4_9_four_planes.dat"));

println("check function1")
function int_gens(I)
    ints = [x for x in I if is_constant(x) && !(x == 0)]
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

for z in 15001:30000
    println(z)
    Mz = matroid_from_revlex_basis_encoding(d4n9[z], 4, 9)
    println("check reference circuits")

    C = [c for c in circuits(Mz) if length(c) == rank(Mz) + 1]

if length(C) == 0
    println("weird")

     open("weird_4_9/4_planes_weird_4_9_2.dat", "a") do file
            write(file,d4n9[z],"\n")
        end

else

    I = matroid_to_reduced_expression_min_chart(Mz, QQ)
    data = join([set_to_stars(I[3],9),d4n9[z]])

 #check realizability   
    if I isa String    
        println("not realizable")
        open("not_realizable_4_9/not_realizable_4_9_2.dat", "a") do file
            write(file,data,"\n")
        end

    elseif isone(ideal(I[1])) 
               println("not realizable unit ideal")

               open("not_realizable_4_9/not_realizable_4_9_2.dat", "a") do file
                     write(file,data,"\n")
                end

    elseif (1 in I[1] || -1 in I[1])
            println("not realizable 1 or -1")
            open("not_realizable_4_9/not_realizable_4_9_2.dat", "a") do file
                write(file,data,"\n")
            end  
#classify if realizable
    elseif (length(I[1]) == 0||I[1] == [0]) #reduces to 0 
            println("zero ideal")
            open("zero_ideal_4_9/reduce_to_zero_4_9_2.dat", "a") do file
            write(file,data,"\n")
        end   
 
    elseif int_gens(I[1])#has integral generator. This only matters if reducing over ZZ

                println("Computing Jacobian ideal with saturation")
                J = simplified_2_singular_locus(I[1],I[2],true)   

            if !(isone(J))
                 println("int gen singular")
                    open("int_generator_4_9/singular_2.dat", "a") do file
                         write(file,data,"\n")
                 end 
            else
                println("int gen smooth")
                    open("int_generator_4_9/smooth_2.dat", "a") do file
                         write(file,data,"\n")
                    end
             end
  

    elseif length(I[1]) == 1
        if (length(vars(I[1][1])) == 1)
                  println("Computing Jacobian ideal with saturation")
                    J = simplified_2_singular_locus(I[1],I[2],true)   
              if !(isone(J))
                     println("princ univariate singular")
                        open("princ_univariate_4_9/singular_2.dat", "a") do file
                             write(file,data,"\n")
                         end 
              else
                    println("princ univariate smooth")
                      open("princ_univariate_4_9/smooth_2.dat", "a") do file
                          write(file,data,"\n")
                        end
              end

            
        else
                 println("Computing Jacobian ideal with saturation")
                 J = simplified_2_singular_locus(I[1],I[2],true)   
             if !(isone(J))
                 println("princ multivariate singular")
                    open("princ_multivariate_4_9/singular_2.dat", "a") do file
                         write(file,data,"\n")
                 end 
            else
                println("princ multivariate smooth")
                    open("princ_multivariate_4_9/smooth_2.dat", "a") do file
                         write(file,data,"\n")
                    end
            end 
        end
    elseif length(I[1])>1
        
        Ivars = ideal_vars(I[1]) 
        if length(Ivars) == 1
                 #println("not principal univariate")#univariate ideal with multiple generators
                        println("Computing Jacobian ideal with saturation")
                      J = simplified_2_singular_locus(I[1],I[2],true)   
                   if !(isone(J))
                     println("not principle univariate singular")
                        open("np_univariate_4_9/singular_2.dat", "a") do file
                             write(file,data,"\n")
                          end 
                    else
                     println("not principle univariate smooth")
                      open("np_univariate_4_9/smooth_2.dat", "a") do file
                          write(file,data,"\n")
                        end
                    end
        else    
            #println("not principle multivariate")#multivariate ideal with multiple generators
                    println("Computing Jacobian ideal with saturation")
                  J = simplified_2_singular_locus(I[1],I[2],true)   
             if !(isone(J))
                 println("np multivariate singular")
                    open("nprinc_multivariate_4_9/singular_2.dat", "a") do file
                         write(file,data,"\n")
                    end 
            else
                println("not principal multivariate smooth")
                    open("nprinc_multivariate_4_9/smooth_2.dat", "a") do file
                         write(file,data,"\n")
                    end
             end
        end
        
    else
       println("to net")#to make sure sorting is exhaustive
       open("net_4_9_2.dat", "a") do file
            write(file,data,"\n")
     end
end
end
end
