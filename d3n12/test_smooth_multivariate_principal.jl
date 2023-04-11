# ARGS[1] = input_file_path
# ARGS[2] = output_file_label




using Combinatorics
using Oscar

currentDir = pwd()

include(joinpath(currentDir, "src/matroid_realization.jl"));
include(joinpath(currentDir, "src/reduceIdealsv2.jl"));
include(joinpath(currentDir, "src/JacobianCriterion.jl"));

d3n12_pm = vec(readlines(joinpath(currentDir, ARGS[1])));
#d3n12_pu = vec(readlines(joinpath(currentDir, "d3n12/sorting_7/principal_univariate.dat")));

function make_directory(dir::AbstractString)
    if !isdir(dir)
        mkdir(dir)
    end
    return dir
end

function simplified_2_singular_locus_with_saturation_check(Igens, Sgens, Qstr, f="")
    L = remove_excess_vars(Igens, Sgens)
    
    I = stepwise_saturation(ideal(L[3]), Sgens)
    
    if length(f) > 0
        if isone(I)  
            open(f, "a")do file
            write(file, Qstr, "\n")
            end
        end
    end
    
    SingL = singular_locus(L[1], L[2], I, L[4]) 
    if SingL isa String
        return SingL
    end
    #J = stepwise_saturation(SingL, L[4])
    return SingL
end

singularDir = make_directory(joinpath(currentDir, "d3n12/singular"))


file_pm_singular = open(joinpath(singularDir, string("principal_multivariate_singular.",ARGS[2],".dat")), "w")
file_pm_smooth = open(joinpath(singularDir, string("principal_multivariate_smooth.",ARGS[2],".dat")), "w")

filename_nonrealizable = joinpath(currentDir, "d3n12",  string("new_nonrealizable.",ARGS[2],".dat"))

filename_skipped_geq3vars = joinpath(singularDir, string("principal_multivariate_geq3vars.",ARGS[2],".dat"))


for z in 1:length(d3n12_pm)

    Mzstr = d3n12_pm[z]
    Mz = matroid_from_revlex_basis_encoding(Mzstr, 3, 12);
    Igens, Sgens = matroid_to_reduced_expression(Mz, QQ);
    
    length(Igens) == 1 || error("not principal")
    
    f = Igens[1]
    vf = vars(f)
    
    if vf > 2
        open(filename_skipped_geq3vars, "a")do file
            write(file, Mzstr, "\n")
        end
    continue
    end
    
    
    I = simplified_2_singular_locus_with_saturation_check(Igens, Sgens, Mzstr, filename_nonrealizable);
    
    
    
    if !isone(I)
        write(file_pm_singular, String(Mstr), "\n")
    else
        write(file_pm_smooth, string(z), "\n")
    end
    println(z, ": ", I)
end

close(file_pm_singular)
close(file_pm_smooth)
