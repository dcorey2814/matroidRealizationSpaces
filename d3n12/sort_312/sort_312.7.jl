
using Combinatorics
using Oscar

currentDir = pwd()

include(joinpath(currentDir, "src/matroid_realization.jl"));
include(joinpath(currentDir, "src/reduceIdealsv2.jl"));


d3n12 = vec(readlines(joinpath(currentDir, "d3n12/3lines_d3n12.7.dat")));


#checks for integer generator

function int_gens(I)
    ints = [x for x in I if is_constant(x)]
    return length(ints)>0    
end

function make_directory(dir::AbstractString)
    if !isdir(dir)
        mkdir(dir)
    end
    return dir    
end 


sortingDir = joinpath(currentDir, "d3n12/sorting_7") # makes the directory "sorting", if not already created.
make_directory(sortingDir)

file_not_realizable = open(joinpath(sortingDir, "not_realizable.dat"), "a") 
file_zero_ideal = open(joinpath(sortingDir, "zero_ideal.dat"), "a")
file_int_generator = open(joinpath(sortingDir, "int_generator.dat"), "a")
file_principal_univariate = open(joinpath(sortingDir, "principal_univariate.dat"), "a")
file_principal_multivariate = open(joinpath(sortingDir, "principal_multivariate.dat"), "a")
file_not_principal_univariate = open(joinpath(sortingDir, "not_principal_univariate.dat"), "a")
file_not_principal_multivariate = open(joinpath(sortingDir, "not_principal_multivariate.dat"), "a")
file_net = open(joinpath(sortingDir, "net.dat"), "a")

#classification over QQ


for z in 1:length(d3n12)

    Mzstr = d3n12[z]
    Mz = matroid_from_revlex_basis_encoding(Mzstr, 3, 12)
    ISdata = matroid_to_reduced_expression(Mz, QQ)
    
    if ISdata isa String    
        write(file_not_realizable, String(Mzstr),"\n"); #println(z, ": not realizable")
        continue
    end
    
    Igens, Sgens = ISdata
    
    if length(Igens) == 0 #reduces to 0 
        write(file_zero_ideal, String(Mzstr), "\n"); #println(z, ": zero")
        continue  
          
    elseif (1 in Igens || -1 in Igens)
        write(file_not_realizable, String(Mzstr),"\n"); #println(z, ": not realizable")
        continue
        
    elseif int_gens(Igens)
        write(file_int_generator, String(Mzstr),"\n"); ; #println(z, ": int_gen")
        continue
        
    elseif length(Igens) == 1
        if (length(vars(Igens[1])) == 1)
            write(file_principal_univariate, String(Mzstr),"\n"); #println(z, ": princ u")
            continue
        else
            write(file_principal_multivariate, String(Mzstr),"\n"); #println(z, ": princ m")
            continue
        end
        
    elseif length(Igens)>1
        
        Ivars = ideal_vars(Igens) 
        if length(Ivars) == 1
            if isone(ideal(Igens)) 
                write(file_not_realizable, String(Mzstr),"\n"); #println(z, ": not realizable")
                continue
            else
                write(file_not_principal_univariate, String(Mzstr),"\n"); #println(z, ": not princ u")
                continue
            end
        else
            write(file_not_principal_multivariate, String(Mzstr),"\n"); #println(z, ": not princ m")
            continue
         end
    else
        write(file_net, String(Mzstr),"\n")
        continue
    end
end

close(file_not_realizable)
close(file_zero_ideal)
close(file_int_generator)
close(file_principal_univariate)
close(file_principal_multivariate)
close(file_not_principal_univariate)
close(file_not_principal_multivariate)
close(file_net)
