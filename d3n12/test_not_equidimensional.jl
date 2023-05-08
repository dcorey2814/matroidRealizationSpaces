# ARGS[1] = input_file_path
# ARGS[2] = output_file_label

using Combinatorics
using Oscar

currentDir = pwd()

include(joinpath(currentDir, "src/matroid_realization.jl"));
include(joinpath(currentDir, "src/reduceIdealsv2.jl"));

d3n12 = vec(readlines(joinpath(currentDir, ARGS[1])));

function make_directory(dir::AbstractString)
    if !isdir(dir)
        mkdir(dir)
    end
    return dir
end

min_prime_dir = make_directory(joinpath(currentDir, "d3n12/min_primes")) 

filename_equidimensional = joinpath(min_prime_dir, string("equidimensional.",ARGS[2],".dat"))
filename_not_equidimensional = joinpath(min_prime_dir, string("not_equidimensional.",ARGS[2],".dat"))


function to_star0(S,n)
    return join(map(x -> x in S ? "*" : "0", 1:n ))
end

io = open(filename_equidimensional, "w")

for z in 1:length(d3n12)

    Mzstr = d3n12[z]
    Mz = matroid_from_revlex_basis_encoding(Mzstr, 3, 12);
    
    charts = maximal_circuits(Mz)
    A = argmax(c -> count_nonbases_chart_int2(Mz, c) , charts)

    Mzstr = to_star0(A,12) * Mzstr
    
    Igens, Sgens, A = matroid_with_chart_to_reduced_expression(Mz, A, QQ);

    length(Igens) == 0 && continue
    any([is_unit(a) for a in Igens]) && continue


    varsIgens = unique!(vcat([vars(f) for f in Igens]...))
    
    I = stepwise_saturation(ideal(Igens), Sgens)
    
    m_primes = minimal_primes(I); 
    not_unique = codim.(m_primes);
    codim_m_primes = unique!(not_unique); 
    println(z,  " numvars : ", length(varsIgens), " gens: ", length(Igens), " codims primes: ", not_unique)
    
    if length(codim_m_primes) > 1
        open(filename_not_equidimensional, "a") do file
            write(file, Mzstr, "\n")
        end
    else
        write(io, Mzstr, "\n")
    end
    
    
    
end

close(io)
