#check all 4-9 matroids s.t original computation gave non principal multivariate ideals. Want to see if there exists a circuit such that the ideal is principal or zero


include("../../src/matroid_realization.jl")
include("../../src/Jacobian_Criterion.jl")

np_m = vec(readlines("not_princ_multivariate_4_9/not_prince_multivariate_4_9.dat"))


#convert circuit to stars and 0s
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

#convert string to reference circuit
function stars_to_set(string,n)

    return Vector{Int64}([a for a in 1:n if string[a] == '*'])
    
end

#computation
for t in 1:length(np_m)
        
   println(t)
    
    Q = np_m[t][10:135]

    Mt = matroid_from_revlex_basis_encoding(Q,4,9)
    
    As = [c for c in circuits(Mt) if length(c) == 5]

    for A in As

        data = join([set_to_stars(A,9),Q])

        MR = new_matroid_realization_space(Mt,A;F = QQ,saturate = true)

        MRR = reduce_ideal_full(MR)

        I = MRR.defining_ideal

        Igens = gens(I)

        if iszero(I)
            println("zero")
            open("recheck/np_m_to_zero.dat", "a") do file
                write(file,data,"\n")
            end
            break

        elseif length(Igens) == 1   
                print("")

                if length(ideal_vars(Igens)) == 1

                    println("princ uni")
                         open("recheck/np_m_to_princ_uni.dat", "a") do file
                            write(file,data,"\n")
                           end
                        break
                else

                       println("pinc mulivar")
                            open("recheck/np_m_to_princ_multi.dat", "a") do file
                            write(file,data,"\n")
                        end
                        break

                end
    end
       end

end

                            


