# creates the jacobian matrix from chosen generators
function jacobian_matrix(R::MPolyRing, x::Vector{<: MPolyElem}, Igens::Vector{<:MPolyElem})
    return matrix(R, [derivative(f,j) for f in Igens, j in 1:length(x) ])
end

# creates the ideal of the singular locus of the ideal I
function singular_locus(R::MPolyRing, x::Vector{<:MPolyElem}, I, Sgens, do_saturation = false) 
    
    if do_saturation
        I = stepwise_saturation(I, Sgens)
    end
        
    m_primes = minimal_primes(I); 
    codim_m_primes = unique!(codim.(m_primes)); 
    
    if length(codim_m_primes) > 1
        return "Not pure codimension"
    end
    c=codim_m_primes[1]
    if c > length(x)
    	return ideal([R(1)])
    end
    Igens = gens(I); 
    Jmat = jacobian_matrix(R, x, Igens); 
    
    Sing = I + ideal(minors(Jmat, c))
    Sing = stepwise_saturation(Sing, Sgens)
    return ideal(R, collect(groebner_basis(Sing)))
end


# after doing the reductions, run this to find the singular locus
function realization_space_2_singular_locus(RS::MatroidRealizationSpace)
    
    !RS.representable && return "not representable"
    
    R= RS.ambient_ring
    x = gens(R)
    I = RS.defining_ideal
    S = RS.inequations
    
    I = stepwise_saturation(I, S)
    isone(I)  && return "not representable"
    
    S = singular_locus(R,x,I,S,true)
end

