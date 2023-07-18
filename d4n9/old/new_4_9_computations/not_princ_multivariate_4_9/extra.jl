using Combinatorics
pm = Polymake

matroids = vec(readlines("not_princ_multivariate_4_9.dat"))

include("../src/matroid_realization.jl")
include("../src/Jacobian_Criterion.jl")

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
println("starting computation")

for m in 1:length(matroids) 

    println(m)

    M = matroid_from_revlex_basis_encoding(matroids[m][10:135],4,9)

    As = optimal_circuits(M)

    has = []

    for A in As

        MR = new_matroid_realization_space(M,A;F = Q,saturate = true)

        MRR = reduce_ideal_full(MR)

        I = MRR.defining_ideal

        if (length(gens(I))>2 || iszero(I) )

                push!(has,A)

        end

    end


    if length(has)>0

            println(has[1])

            data = join([has[1],matroids[m][10:135]])


            variables = ideal_vars(Igens)


    

                
