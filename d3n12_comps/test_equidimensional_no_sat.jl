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

filename_zero = joinpath(min_prime_dir, string("zero.",ARGS[2],".dat"))
filename_one = joinpath(min_prime_dir, string("one.",ARGS[2],".dat"))
filename_univariate = joinpath(min_prime_dir, string("univariate.",ARGS[2],".dat"))
filename_principal = joinpath(min_prime_dir, string("principal.",ARGS[2],".dat"))

filename_equidimensional = joinpath(min_prime_dir, string("var2_equidimensional.",ARGS[2],".dat"))
filename_not_equidimensional = joinpath(min_prime_dir, string("var2_not_equidimensional.",ARGS[2],".dat"))

filename_var3 = joinpath(min_prime_dir, string("vargeq3.",ARGS[2],".dat"))

function to_star0(S,n)
    return join(map(x -> x in S ? "*" : "0", 1:n ))
end


io_zero = open(filename_zero, "a") 
io_one = open(filename_one, "a") 
io_uni = open(filename_univariate, "a") 
io_prin = open(filename_principal, "a")
io_var2 = open(filename_equidimensional, "a")
io_var3 = open(filename_var3, "a")


function matroid_with_chart_to_reduced_expression_no_saturation(Q, A, F)
     
    RQ = matroid_realization_space(Q, A, F)
    R = parent(RQ[1][1])
    Sgens = RQ[2]    
    Sgens = gens_2_factors(Sgens)
    
    Igens_notsat = gens(ideal(RQ[1]))
    
    reducedData = reduce_ideal_full(Igens_notsat, Sgens, R, gens(R), false)
    reducedData isa String && return reducedData
    
    Igens = reducedData[1]
    Sgens = reducedData[2]
    
    if length(Igens) == 0 
        Igens = [R(0)]
    end 

    Igens = collect(groebner_basis(ideal(Igens)))
        
    return (Igens, Sgens, A)
    
end


function small_var3(Igens)
    varg


end

for z in parse(Int64, ARGS[3]):length(d3n12)

    Mzstr = d3n12[z]
    Mz = matroid_from_revlex_basis_encoding(Mzstr, 3, 12);
    
    charts = maximal_circuits(Mz)
    mx = maximum(c -> count_nonbases_chart_int2(Mz, c) , charts)
    options = [c for c in charts if count_nonbases_chart_int2(Mz, c) >= mx-1]
    
    if length(options) >=7
        A = options[7]
    elseif length(charts) > 6
                println("not enough max charts")
        A = charts[7]
    else
        println("not enough charts")
        A = last(charts)
    end
        
    
    Mzstr = to_star0(A,12) * Mzstr
    
    Igens, Sgens, A = matroid_with_chart_to_reduced_expression_no_saturation(Mz, A, QQ);


    if length(Igens) == 0 
        write(io_zero, Mzstr, "\n")
        continue
    elseif any([is_unit(a) for a in Igens])
        write(io_one, Mzstr, "\n")
        continue
    end
    
    
    varsIgens = unique!(vcat([vars(f) for f in Igens]...))
    
    println(z,  " numvars : ", length(varsIgens), " gens: ", length(Igens))
    
    if length(Igens) == 0 
        write(io_zero, Mzstr, "\n")
        continue
    elseif any([is_unit(a) for a in Igens])
        write(io_one, Mzstr, "\n")    
        continue
    elseif length(varsIgens) == 1
        write(io_uni, Mzstr, "\n")
        continue
    elseif length(Igens) == 1
        write(io_prin, Mzstr, "\n")
        continue
    end
    
    if length(varsIgens) <= 2
        I = stepwise_saturation(ideal(Igens), Sgens)
        Igens = gens(I)
        
        m_primes = minimal_primes(I); 
        codim_m_primes = codim.(m_primes);
        println(z, " codims primes: ", codim_m_primes)
        unique!(codim_m_primes); 
    
        if length(codim_m_primes) > 1
            open(filename_not_equidimensional, "a") do file
                write(file, Mzstr, "\n")
            end
            continue
    
        else
            write(io_var2, Mzstr, "\n")
        end
    else
        write(io_var3, Mzstr, "\n")
    end    
end

close(io_zero)
close(io_one)
close(io_uni)
close(io_prin)
close(io_var2)
close(io_var3)