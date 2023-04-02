using Combinatorics
using Oscar

currentDir = pwd()


include(joinpath(currentDir, "src/matroid_realization.jl"));
include(joinpath(currentDir, "src/reduceIdealsv2.jl"));


d3n10_3_lines = vec(readlines(joinpath(currentDir, "d3n10/3lines_d3n10.dat")));

#checks for integer generator

function int_gens(I)
    ints = [x for x in I if is_constant(x)]
    return length(ints)>1    
end

function make_directory(dir::AbstractString)
    if !isdir(dir)
        mkdir(dir)
    end
    return dir    
end 


sortingDir = joinpath(currentDir, "d3n10/sorting")
make_directory(sortingDir)

file_not_realizable = open(joinpath(sortingDir, "not_realizable.dat"), "w") 
file_zero_ideal = open(joinpath(sortingDir, "zero_ideal.dat"), "w")
file_int_generator = open(joinpath(sortingDir, "int_generator.dat"), "w")
file_principal_univariate = open(joinpath(sortingDir, "principal_univariate.dat"), "w")
file_principal_multivariate = open(joinpath(sortingDir, "principal_multivariate.dat"), "w")
file_not_principal = open(joinpath(sortingDir, "not_principal.dat"), "w")
file_net = open(joinpath(sortingDir, "net.dat"), "w")

#classification over QQ

for z in 1:length(d3n10_3_lines)

    Mzstr = d3n10_3_lines[z]
    Mz = matroid_from_revlex_basis_encoding(Mzstr, 3, 10)
    ISdata = matroid_to_reduced_expression(Mz, QQ)
    
    if ISdata isa String    
        write(file_not_realizable, String(Mzstr),"\n")
        continue
    end
    
    Igens, Sgens = ISdata
    
    if length(Igens) == 0 #reduces to 0 
        write(file_zero_ideal, String(Mzstr), "\n")
        continue  
          
    elseif (1 in Igens || -1 in Igens)
        write(file_not_realizable, String(Mzstr),"\n")
        continue
        
    elseif int_gens(Igens)
        write(file_int_generator, String(Mzstr),"\n")
        continue
        
    elseif length(Igens) == 1
        if (length(vars(Igens[1])) == 1)
            write(file_principal_univariate, String(Mzstr),"\n")
            continue
        else
            write(file_principal_univariate, String(Mzstr),"\n")
            continue
        end
        
    elseif length(Igens)>1
        write(file_not_principal, String(Mzstr),"\n")
        continue
        
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
close(file_not_principal)
close(file_net)
