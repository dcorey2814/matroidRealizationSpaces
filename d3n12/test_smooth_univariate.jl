using Combinatorics
using Oscar

currentDir = pwd()

include(joinpath(currentDir, "src/matroid_realization.jl"));
include(joinpath(currentDir, "src/reduceIdealsv2.jl"));
include(joinpath(currentDir, "src/JacobianCriterion.jl"));


d3n12_pu = vec(readlines(joinpath(currentDir, "d3n12/sorting_7/principal_univariate.dat")));


#checks for integer generator


function make_directory(dir::AbstractString)
    if !isdir(dir)
        mkdir(dir)
    end
    return dir
end

singularDir = make_directory(joinpath(currentDir, "d3n12/singular"))

file_pm_singular = open(joinpath(singularDir, "principal_univariate_singular.dat"), "w")
file_pm_smooth = open(joinpath(singularDir, "principal_univariate_smooth.dat"), "w")


for z in 1:length(d3n12_pu)

    Mzstr = d3n12_pu[z]
    Mz = matroid_from_revlex_basis_encoding(Mzstr, 3, 12)
    Igens, Sgens = matroid_to_reduced_expression(Mz, QQ)

    I = simplified_2_singular_locus(Igens, Sgens)
    
    if !isone(I)
        write(file_pm_singular, String(Mstr), "\n")
    else
        write(file_pm_smooth, string(z), "\n")
    end
    println(z, ": ", I)

    #I = simplified_2_singular_locus(Igens, Sgens)
    #if !isone(I)
    #    write(file_pm_singular, String(Mzstr), "\n")
    #else
    #    write(file_pm_smooth, String(Mzstr), "\n")
    #end
    #println(z, ": ", length(Igens[1]))

end

close(file_pm_singular)

