using Oscar
using Combinatorics
pm = Polymake

currentDir = pwd()

d4n9 = vec(readlines("../data/d4n9.dat"));

include("../src/matroid_realization.jl");
include("../src/reduceIdealsv2.jl")
include("../src/JacobianCriterion.jl")

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

#check not principal multivariate
npmult = vec(readlines("sorting/not_principal_multivariate.dat"))
#io1npmult = open(joinpath(pwd(), "test_nprinc_mult_singular_4_9.dat"), "w")
#io2npmult = open(joinpath(pwd(), "test_nprinc_mult_smooth_4_9.dat"), "w")

for t in 7618:8000 

    println(t)
    Mt = matroid_from_revlex_basis_encoding(npmult[t],4,9)

    println("reducing ideal")
    I = matroid_to_reduced_expression_fewest_var(Mt,QQ)
    data = join([set_to_stars(I[3],9),npmult[t]])
    
    if length(ideal_vars(I[1]))>3
        println("recheck with different circuit")
        open("to_recheck_nprinc_mult_singular_4_9.dat", "a") do file
                write(file,data,"\n")
            end
    else

        println("checking jacobian criteria")
        J = simplified_2_singular_locus(I[1], I[2])
            if !(is_one(J))
                println("                  singular")
                open("final_nprinc_mult_singular_4_9.dat", "a") do file
                    write(file,data,"\n")
                end
            else
                println("                  smooth")
                open("final_nprinc_mult_smooth_4_9.dat", "a") do file
                     write(file,data,"\n")
                end
            end

         end
       end


#close(io1npmult)
#close(io2npmult)
#remember 7610 got checke twice, find way to get duplicate out of file
