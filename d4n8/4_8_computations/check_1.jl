using Oscar
using Combinatorics
pm = Polymake

#stuck on: 52 (jacobian ideal) #105
#8270 stuck on realization space

include("../../src/matroid_realization.jl")
include("../../src/Jacobian_Criterion.jl")

#load 4-8 matroids satisfying 4-lines property
d4n8 = vec(readlines("4_8_four_planes_connected.dat"));



#converts set to stars and 0s string
println("check function 1")
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


#find reference circuit that gives simplest ideal. This seemed to work better for interacting with the jacobian code
function count_nonbases_disjoint_to_chart(Q, A)
    NBs = nonbases(Q)
    return length([nb for nb in NBs if length(intersect(A,nb)) == 0])
end



println("start loop")
for z in 1:length(d4n8)
    println(z)
    Mz = matroid_from_revlex_basis_encoding(d4n8[z], 4, 8)
    
    println("check connectivity")

    if !is_connected(Mz)

        println("disconnected")

        open("disconnected_4_8_1.dat", "a") do file
            write(file,d4n8[z],"\n")
          end
    else
    println("check reference circuits")#check if reference circuit exists

    C = [c for c in circuits(Mz) if length(c) == rank(Mz) + 1]

    if length(C) == 0
        println("no ref circ")

         open("4_8_no_ref_circ.dat", "a") do file
            write(file,d4n8[z],"\n")
          end

    else #contintinue if reference circuit exists
    println("compute realization space")   

    #find optimal reference circuit

    #As = optimal_circuits(Mz)
    
    A = argmin(c -> count_nonbases_disjoint_to_chart(Mz, c) , C)

    println(A)

    data = join([set_to_stars(A,8),d4n8[z]])#data to be recorded

    MR = new_matroid_realization_space(Mz,A;F = QQ,saturate = true) #compute realization space


        if !MR.representable#check realizability

            println("not realizable")
            open("4_8_not_realizable.dat", "a") do file
                write(file,data,"\n")
            end
           
        else  #if realizable, reduce data

                println("reducing ring")
                MRR = reduce_ideal_full(MR)
    
                I = MRR.defining_ideal

                Igens = gens(I)
          
        #classify reduced data
            println("classify")
            if iszero(I)#check if ideal reduces to zero wrt to this reference circuit
            
                println("zero ideal")
                    open("4_8_zero.dat", "a") do file
                    write(file,data,"\n")
            
                end 
    

            elseif length(Igens) == 1 #check if principal

                    if length(ideal_vars(Igens)) == 1#check if univariate
                        println("computing jacobian ideal")
                        J = realization_space_2_singular_locus(MRR)#check if singular

                        if isone(J)#smooth
                          println("princ univariate smooth")
                          open("4_8_princ_univariate_smooth.dat", "a") do file
                             write(file,data,"\n")
                         end 
                    else #singular

                        println("princ univariate singular")
                        open("4_8_princ_univariate_singular.dat", "a") do file
                          write(file,data,"\n")
                        end
                    end

                 else   #then its multivariate
                       println("computing jacobian ideal") 
                    J = realization_space_2_singular_locus(MRR)#check if singular

                    if isone(J)#smooth
                           println("princ multivariate smooth")
                          open("4_8_princ_multivariate_smooth.dat", "a") do file
                             write(file,data,"\n")
                         end 
                    else #singular

                        println("princ multivariate singular")
                        open("4_8_princ_multivariate_singular.dat", "a") do file
                          write(file,data,"\n")
                        end


                    end


                end        
         elseif length(Igens) > 1 #check if multiple generators

                    if length(ideal_vars(Igens)) == 1#check if univariate
                        
                        println("computing jacobian ideal")

                        J = realization_space_2_singular_locus(MRR)#check if singular

                      if isone(J)#smooth
                              println("not princ univariate smooth")
                              open("4_8_not_princ_univariate_smooth.dat", "a") do file
                                 write(file,data,"\n")
                             end 
                      else #singular
    
                            println("not princ univariate singular")
                            open("4_8_not_princ_univariate_singular.dat", "a") do file
                              write(file,data,"\n")
                            end

                        end

    
                    else#then its multivariate
                        
                        println("computing jacobian ideal")
                        J = realization_space_2_singular_locus(MRR)#check if singular

                        if isone(J)#smooth
                               println("not princ multivariate smooth")
                              open("4_8_not_princ_multivariate_smooth.dat", "a") do file
                                 write(file,data,"\n")
                             end 
                        else #singular

                            println("not multivariate singular")
                            open("4_8_not_princ_multivariate_singular.dat", "a") do file
                              write(file,data,"\n")
                            end

                        end

                   end

             else
       
                println("to net")#to make sure sorting is exhaustive
                open("net_4_8.dat", "a") do file
                 write(file,data,"\n")
                end
         end
end
end
end
end
                    


           

        





































  
