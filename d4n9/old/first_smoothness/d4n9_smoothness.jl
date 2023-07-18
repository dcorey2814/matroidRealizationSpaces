using Oscar
using Combinatorics
pm = Polymake

currentDir = pwd()

d4n9 = vec(readlines("../data/d4n9.dat"));

include("../src/matroid_realization.jl");
include("../src/reduceIdealsv2.jl")
include("../src/JacobianCriterion.jl")
include("../src/TSC_CoordRingV2.jl")


sp_49 = load("4planes_4_9.dat")


pmult = vec(readlines("sorting/principal_multivariate.dat"))
puni = vec(readlines("sorting/principal_univariate.dat"))
npmult = vec(readlines("sorting/not_principal_multivariate.dat"))
npuni = vec(readlines("sorting/not_principal_univariate.dat"))  


#check principal multivariate

pmult = vec(readlines("sorting/principal_multivariate.dat"))
io1pmult = open(joinpath(pwd(), "2_princ_mult_singular_4_9.dat"), "w")
io2pmult = open(joinpath(pwd(), "2_princ_mult_smooth_4_9.dat"), "w")

for t in 1:length(pmult)
    Mt = matroid_from_revlex_basis_encoding(pmult[t],4,9)
    I = matroid_to_reduced_expression(Mt,QQ)
   
    J = simplified_2_singular_locus(I[1], I[2])
        if !(is_one(J))
            open(io1pmult, "a") do file
                write(file, String(pmult[t]),"\n")
            end
        else
            open(io2pmult, "a") do file
                 write(file, String(pmult[t]),"\n")
            end
        end
        
end
close(io1pmult)
close(io2pmult)



#check principal univariate
puni = vec(readlines("sorting/principal_univariate.dat"))
io1puni = open(joinpath(pwd(), "2_princ_uni_singular_4_9.dat"), "w")
io2puni = open(joinpath(pwd(), "2_princ_uni_smooth_4_9.dat"), "w")

for t in 1:length(puni)
    Mt = matroid_from_revlex_basis_encoding(puni[t],4,9)
    I = matroid_to_reduced_expression(Mt,QQ)
   
    J = simplified_2_singular_locus(I[1], I[2])
        if !(is_one(J))
            open(io1puni, "a") do file
                write(file, String(puni[t]),"\n")
            end
        
        else
        
            open(io1puni, "a") do file
                write(io2puni, String(puni[t]),"\n")
            end
         end
end
close(io1puni)
close(io2puni)


#check not principal multivariate
npmult = vec(readlines("sorting/not_principal_multivariate.dat"))
io1npmult = open(joinpath(pwd(), "nprinc_mult_singular_4_9.dat"), "w")
io2npmult = open(joinpath(pwd(), "nprinc_mult_smooth_4_9.dat"), "w")

for t in 1:length(npmult)
    Mt = matroid_from_revlex_basis_encoding(npmult[t],4,9)
    I = matroid_to_reduced_expression(Mt,QQ)
   
    J = simplified_2_singular_locus(I[1], I[2])
        if !(is_one(J))
            open(io1pnmult, "a") do file
                write(file, String(npmult[t]),"\n")
            end
        else
        
            open(io2pnmult, "a") do file
                 write(file, String(npmult[t]),"\n")
            end
        end
   end

close(io1npmult)
close(io2npmult)

#check not principal univariate
npuni = vec(readlines("sorting/not_principal_univariate.dat"))  
io1npuni = open(joinpath(pwd(), "2_nprinc_uni_singular_4_9.dat"), "w")
io2npuni = open(joinpath(pwd(), "2_nprinc_uni_smooth_4_9.dat"), "w")

for t in 1:length(npuni)
    Mt = matroid_from_revlex_basis_encoding(npuni[t],4,9)
    I = matroid_to_reduced_expression(Mt,QQ)
   
    J = simplified_2_singular_locus(I[1], I[2])
        if !(is_one(J))
        
            open(io1npuni, "a") do file
                write(file, String(npuni[t]),"\n")
            end
        else
             open(io2npuni, "a") do file
                 write(file, String(npuni[t]),"\n")
            end
         end
        
end
close(io1npuni)
close(io2npuni)


