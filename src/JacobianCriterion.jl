
# creates the jacobian matrix from chosen generators
function jacobian_matrix(R::MPolyRing, x::Vector{T}, Igens::Vector{T}) where T <: MPolyElem
    return matrix(R, [derivative(f,j) for f in Igens, j in 1:length(x) ])
end


# creates the ideal of the singular locus of the ideal I
function singular_locus(R::MPolyRing, x::Vector{T}, I) where T <: MPolyElem
    m_primes = minimal_primes(I); 
    codim_m_primes = unique!(codim.(m_primes)); 
    if length(codim_m_primes) > 1
        error("Not pure codimension")
    end
    c=codim_m_primes[1]
    if c > length(x)
    	return ideal([R(1)])
    end
    Igens = gens(I); 
    Jmat = jacobian_matrix(R, x, Igens) ;
    return I + ideal(minors(Jmat, c))
end



# realizes ideal and semigroup generators in a smaller polynomial ring
function remove_excess_vars(Igens, Sgens)
    vars_used = union!([vars(f) for f in union(Igens, Sgens)]...)
    
    R = parent(Igens[1])
    varlist = gens(R)
    Rred, new_varlist = PolynomialRing(coefficient_ring(R) , :y=>(1:length(vars_used)))
    pr = []; 
    j=1
    for v in varlist
        if v in vars_used
            push!(pr, new_varlist[j])
            j=j+1
        else
            push!(pr, Rred(0))
        end
    end
    
    phi = hom(R, Rred, a->a, pr)
    return (Rred, new_varlist, phi.(Igens), phi.(Sgens))
end


# after doing the reductions, run this to find the singular locus
function simplified_2_singular_locus(Igens, Sgens)
    L = remove_excess_vars(Igens, Sgens)
    SingL = singular_locus(L[1], L[2], ideal(L[3])) 
    S = MPolyPowersOfElement(L[1] , L[4])
    SinvR, inc = Localization(S)
    return saturated_ideal(inc(SingL))
end
