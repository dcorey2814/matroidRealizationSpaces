using Oscar
using Combinatorics
pm = Polymake

currentDir = pwd()

d4n9 = vec(readlines("../data/d4n9.dat"));

include("../src/matroid_realization.jl");
include("../src/reduceIdealsv2.jl")
include("../src/JacobianCriterion.jl")

function set_to_stars(set)
    
    z = []
    
    for a in 1:12
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

for t in 1:100
    Mt = matroid_from_revlex_basis_encoding(npmult[t],4,9)
    #I = matroid_to_reduced_expression_fewest_var(Mt,QQ)
    I = matroid_to_reduced_expression(Mt,QQ)

    #println(" ",t," with ", I[3])
    println(t)
    #data = join([set_to_stars(I[3]),npmult[t]])
    data = npmult[t]
   
    J = simplified_2_singular_locus(I[1], I[2])
        if !(is_one(J))
            println("                  singular")
            open("test_59_nprinc_mult_singular_4_9.dat", "a") do file
                write(file,data,"\n")
            end
        else
            println("                  smooth")
            open("test_59_nprinc_mult_smooth_4_9.dat", "a") do file
                 write(file,data,"\n")
            end
        end
   end

#close(io1npmult)
#close(io2npmult)
