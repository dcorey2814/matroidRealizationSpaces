using Oscar
using Combinatorics
pm = Polymake

#classify 3-12 matroids with singular realization spaces
include("../src/matroid_realization.jl")
include("../src/Jacobian_Criterion.jl")

#load singular 3-12
matroids = vec(readlines("singular_3_12.dat"));



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
for z in 7:76
    println(z)
    Mz = matroid_from_revlex_basis_encoding(matroids[z], 3, 12)
    
    println("check connectivity")

    if !is_connected(Mz)

        println("disconnected")

        open("disconnected_3_12/disconnected_3_12_1.dat", "a") do file
            write(file,matroids[z],"\n")
          end
    else
    println("check reference circuits")#check if reference circuit exists

    C = [c for c in circuits(Mz) if length(c) == rank(Mz) + 1]

    if length(C) == 0
        println("no ref circ")

         open("no_ref_circ_3_12/no_ref_circ_1.dat", "a") do file
            write(file,matroids[z],"\n")
          end

    else #contintinue if reference circuit exists
    println("compute realization space")   

    #find optimal reference circuit

    #As = optimal_circuits(Mz)
    
    A = optimal_circuits(Mz)[1]

    println(A)

    data = join([set_to_stars(A,12),matroids[z]])#data to be recorded

    MR = new_matroid_realization_space(Mz,A;F = QQ,saturate = true) #compute realization space


        if !MR.representable#check realizability

            println("not realizable")
            open("not_realizable_3_12_not_realizable_3_12_1.dat", "a") do file
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
                    open("zero_test.dat", "a") do file
                    write(file,data,"\n")
            
                end 
    

            elseif length(Igens) == 1 #check if principle

                    if length(ideal_vars(Igens)) == 1#check if univariate
                        println("computing jacobian ideal")
                        J = realization_space_2_singular_locus(MRR)#check if singular

                        if isone(J)#smooth
                          println("princ univariate smooth")
                          open("princ_univariate_3_12_smooth_1.dat", "a") do file
                             write(file,data,"\n")
                         end 
                    else #singular

                        println("princ univariate singular")
                        open("princ_univariate_3_12_singular_1.dat", "a") do file
                          write(file,data,"\n")
                        end
                    end

                 else   #then its multivariate
                       println("computing jacobian ideal") 
                    J = realization_space_2_singular_locus(MRR)#check if singular

                    if isone(J)#smooth
                           println("princ multivariate smooth")
                          open("princ_multivariate_3_12_smooth_1.dat", "a") do file
                             write(file,data,"\n")
                         end 
                    else #singular

                        println("princ multivariate singular")
                        open("princ_multivariate_3_12_singular_1.dat", "a") do file
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
                              open("not_princ_univariate_3_12_smooth_1.dat", "a") do file
                                 write(file,data,"\n")
                             end 
                      else #singular
    
                            println("not princ univariate singular")
                            open("not_princ_univariate_3_12_singular_1.dat", "a") do file
                              write(file,data,"\n")
                            end

                        end

    
                    else#then its multivariate
                        
                        println("computing jacobian ideal")
                        J = realization_space_2_singular_locus(MRR)#check if singular

                        if isone(J)#smooth
                               println("not princ multivariate smooth")
                              open("not_princ_multivariate_3_12_smooth_1.dat", "a") do file
                                 write(file,data,"\n")
                             end 
                        else #singular

                            println("not princ multivariate singular")
                            open("not_princ_multivariate_3_12_singular_1.dat", "a") do file
                              write(file,data,"\n")
                            end

                        end

                   end

             else
       
                println("to net")#to make sure sorting is exhaustive
                open("net_3_12_1.dat", "a") do file
                 write(file,data,"\n")
                end
         end
end
end
end
end
                    


           

        





































  
