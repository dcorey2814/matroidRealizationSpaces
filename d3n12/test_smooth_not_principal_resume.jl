# ARGS[1] = input_file_path
# ARGS[2] = output_file_label

using Combinatorics
using Oscar

currentDir = pwd()

include(joinpath(currentDir, "src/matroid_realization.jl"));
include(joinpath(currentDir, "src/reduceIdealsv2.jl"));
include(joinpath(currentDir, "src/JacobianCriterion.jl"));

d3n12 = vec(readlines(joinpath(currentDir, ARGS[1])));

function make_directory(dir::AbstractString)
    if !isdir(dir)
        mkdir(dir)
    end
    return dir
end

function simplified_2_singular_locus_with_saturation_check(Igens, Sgens, Qstr, f="")


    if length(Igens) == 0
        return ideal([1])
    end
    L = remove_excess_vars(Igens, Sgens)

    I = stepwise_saturation(ideal(L[3]), L[4])

    if length(f) > 0
        if isone(I)
            open(f, "a")do file
                write(file, Qstr, "\n")
            end
        return "is one"
        end
    end

    SingL = singular_locus(L[1], L[2], I, L[4])
    if SingL isa String
        return SingL
    end
    return SingL
end

singularDir = make_directory(joinpath(currentDir, "d3n12/singular_not_principal"))

filename_zeroideal = joinpath(singularDir, string("zero_ideal.",ARGS[2],".dat"))

filename_principal_smooth = joinpath(singularDir, string("principal_smooth.",ARGS[2],".dat"))
filename_principal_singular = joinpath(singularDir, string("principal_singular.",ARGS[2],".dat"))


filename_np_2vars_smooth = joinpath(singularDir, string("not_principal_2var_smooth.",ARGS[2],".dat"))
filename_np_2vars_singular = joinpath(singularDir, string("not_principal_2var_singular.",ARGS[2],".dat"))

filename_np_3vars_smooth = joinpath(singularDir, string("not_principal_3var_smooth.10v.", ARGS[2],".dat"))
filename_np_3vars_singular = joinpath(singularDir, string("not_principal_3var_singular.10v.", ARGS[2],".dat"))

filename_np_3vars_remaining = joinpath(singularDir, string("not_principal_3var_remaining.", ARGS[2],".dat"))
filename_np_4vars_remaining = joinpath(singularDir, string("not_principal_geq4var_remaining.", ARGS[2],".dat"))

filename_nonrealizable = joinpath(singularDir, string("new_nonrealizable.npm.",ARGS[2],".dat"))


file_2_sm = open(filename_np_2vars_smooth, "a")
file_3_sm = open(filename_np_3vars_smooth, "a")

file_3 = open(filename_np_3vars_remaining, "a")
file_4 = open(filename_np_4vars_remaining, "a")

for z in parse(Int64, ARGS[3]):length(d3n12)

    Mzstr = d3n12[z]
    Mz = matroid_from_revlex_basis_encoding(Mzstr, 3, 12);
    Igens, Sgens = matroid_to_reduced_expression(Mz, QQ);

    varsIgens = unique!(vcat([vars(f) for f in Igens]...))
    println(z, " numvars : ", length(varsIgens), " gens: ", length(Igens))


    
    if length(Igens) == 0
        open(filename_zeroideal, "a") do file
            write(file, Mzstr, "\n")
        end
        continue
    end 
    
    R = parent(Igens[1])
    
    if Igens == [R(0)] 
        open(filename_zeroideal, "a") do file
            write(file, Mzstr, "\n")
        end
    
    elseif length(Igens) == 1
        I = simplified_2_singular_locus_with_saturation_check(Igens, Sgens, Mzstr, filename_nonrealizable);
        
        if I isa String
            continue
        
        elseif !isone(I)
            open(filename_principal_singular, "a") do file
                write(file, Mzstr, "\n")
            end
        else
            open(filename_principal_smooth, "a") do file
                write(file, Mzstr, "\n")
            end
        end

    elseif length(varsIgens) == 2
        I = simplified_2_singular_locus_with_saturation_check(Igens, Sgens, Mzstr, filename_nonrealizable);
        
        if I isa String
            continue
        end
        
        if !isone(I)
            open(filename_np_2vars_singular, "a") do file
                write(file, Mzstr, "\n")
            end
        else
            write(file_2_sm, Mzstr, "\n")
        end

    elseif length(varsIgens) == 3

        if maximum([length(f) for f in Igens]) <= 10

            I = simplified_2_singular_locus_with_saturation_check(Igens, Sgens, Mzstr, filename_nonrealizable);
            
            if I isa String
                continue
            end
            
            if !isone(I)
                open(filename_np_3vars_singular, "a") do file
                    write(file, Mzstr, "\n")
                end
            else
                write(file_3_sm, Mzstr, "\n")
            end
        else
	    J = stepwise_saturation(ideal(Igens), Sgens)
            if isone(J)
                open(filename_nonrealizable, "a")do file
                    write(file, Mzstr, "\n")
                end
            else
                write(file_3, Mzstr, "\n")
            end
        end
    else
        J = stepwise_saturation(ideal(Igens), Sgens)
        if isone(J)
            open(filename_nonrealizable, "a")do file
                write(file, Mzstr, "\n")
            end
        else
            write(file_4, Mzstr, "\n")
        end
    end
end

close(file_2_sm)
close(file_3_sm)
close(file_3)
close(file_4)
