#code that computes the coordinate ring of a TSC with respect to some matroid, basis, and field. We can also get the associated matrix with respect to a basis and field. Also, we can get the reduced presentation of the coordinate ring given some basis and field.

#Note that this code needs to interact with the code from matroid_realization.jl

#Throughout Q = matroid, A = basis, F= field.

################################################
#functions to reduce ideals and semigrouü


#compute coefficient seperately

function coefficient_v(v,f,R)
    
    (degree(f,v) == 1) || return "degree of variable must be 1"

    withv = [term(f,i) for i in 1:length(f) if v in vars(monomial(f,i))]
    g = sum(withv)
    
    gFactor = factor(g); 
    ug = unit(gFactor); 
    gFactorDict = Dict([pair for pair in gFactor]) ; 
    
    factors_not_v = [k^(gFactorDict[k]) for k in keys(gFactorDict) if k ≠ v ]
    
    if length(factors_not_v) == 0
        return ug
    else
        return ug*prod(factors_not_v)
    end    
end

# returns the factors of f, but not the exponents
function poly_2_factors(f)
    return collect(keys(Dict(factor(f))))
end

# returns the unique factors of the elements of Sgen, again no exponents. 
function gens_2_factors(Sgens)
    return unique!(vcat([poly_2_factors(f) for f in Sgens]...))
end


#function that finds generators of ideal that contain x, checks if coefficient is in semigroup, then solves
function find_solution_v(v, Igens, R, Sall)
    
    with_v_deg_1 = [gen for gen in Igens if (v in vars(gen) && degree(gen,v)==1)] 
    (length(with_v_deg_1) != 0) || return "can't isolate"

    for f in with_v_deg_1
       #println("v: ", v)        
       #println("f: ", factor(f))

        no_v = [term(f,i) for i in 1:length(f) if !(v in vars(monomial(f,i)))]

        if length(no_v) == 0
            continue
        end

#        (length(no_v) != 0) || return "can't isolate, variable is a factor"
        
        den = coefficient_v(v, f, R)
        fac_den = poly_2_factors(den)
       #println("den: ", den)
        
        #issubset(fac_den, Sall) || return "can't isolate, coefficient of variable is not a unit"
        if !issubset(fac_den, Sall)
            continue
        end
        
        h = R(-1)*sum(no_v)
        return h//den    
    
    end
    return "can't solve for v"
end

function sub_map(v, t, R, varlist) # v is replaced by t in f
    sub_x = []
    for i in 1:length(varlist)
        if v == varlist[i]
            push!(sub_x, t)
        else
            push!(sub_x, varlist[i])
        end
    end
    return hom(R,FractionField(R),a->a,sub_x)
end


function sub_v(v, t, f, R, varlist) # replace v by t in f, only return the numerator.
    m = sub_map(v,t,R,varlist) 
    new_f = numerator(m(f))
    return new_f #m(f) #new_f
end

#variables in ideal
function ideal_vars(Igens) 
    return unique!(vcat([vars(gen) for gen in Igens]...))
end

#function to produce new ideal generators
function n_new_Igens(x, tx, Igens, newSgens, R, varlist)    
    preIgens = unique!([clean(sub_v(x, tx, gen, R, varlist), R, newSgens) for gen in Igens])
    return filter(x-> x!= R(0), preIgens)
end

#function to produce new subgroup generators
function n_new_Sgens(x, tx, Sgens, R, varlist)
    preSgens = unique!([sub_v(x, tx, gen, R, varlist) for gen in Sgens])
    return gens_2_factors(preSgens)
end

function reduce_ideal_one_step(Igens, Sgens, R, varlist, fullyReduced)
    
    if R(0) in Sgens
        
        return "Not realizable 0 in Semigroup"
        
    else
        Ivars = ideal_vars(Igens); 
        
        for x in Ivars 
           tx = find_solution_v(x, Igens, R, Sgens)
            if tx isa String
                continue
            else 
                Sgens_new = n_new_Sgens(x, tx, Sgens, R, varlist);
                Igens_new = n_new_Igens(x, tx, Igens, Sgens_new, R, varlist); 

                return (Igens_new, Sgens_new, R, varlist, fullyReduced)
         
            end
        

        end
    return (Igens, Sgens, R, varlist, true)
    end 
end


function reduce_ideal_full(Igens, Sgens, R, varlist, fullyReduced = false)
    
    output = reduce_ideal_one_step(Igens, Sgens, R, varlist, fullyReduced)
    if output isa String
    #if reduce_ideal_one_step(Igens, Sgens, R, varlist, fullyReduced) isa String
        
        return "Not Realizable 0 in Semigroup"
        
    else
    (Igens, Sgens, R, varlist, fullyReduced) = output
        if !fullyReduced
            
            return reduce_ideal_full(Igens, Sgens, R, varlist, fullyReduced)
        else
            return (Igens, Sgens, R, varlist, fullyReduced)
        end
    end
end

#functions used in TSC computation

# given the bases Bs of a matroid, and a fixed basis B, this function finds
# the nonzero coordinates xij of the coordinate ring of the matroid stratum,
# These correspond to all elements A of Bs such that the symmetric difference
# with B has exactly 2 elements. 

function bases_matrix_coordinates(Bs::Vector{Vector{Int}}, B::Vector{Int})
    
    coord_bases = [b for b in Bs if length(symdiff(B,b)) == 2]
    
    new_coords = Vector{Vector{Int}}([])
    
    for b in coord_bases
        row_b = setdiff(B,b)[1]
        row_b = count(a->(a<row_b), B) + 1
        # count(f, v) does what?
        # - iterate through the elements a in v
        # - compute f(a) 
        # - if that is true, increment the counting variable by 1
        # - otherwise, continue
        # - return the value of the internal counter.
        # Similar with all(a->(a<row_b), B), for instance.
        #row_b = length([a for a in B if a < row_b]) + 1
        
        col_b = setdiff(b,B)[1]
        col_b = col_b - length([a for a in B if a ≤ col_b]) 
        
        push!(new_coords, [row_b,col_b])
    end
    return sort!(new_coords, by = x -> (x[2], x[1]))    
end



# Given the bases Bs of a matroid, a fixed basis B, and a coefficient field F
# this function creates a polynomial ring in xij, where the xij are 
# determined by the function basis_matrix_coordinates. This function also
# returns (as the 3rd element of a triple) a dictionary (i,j) => xij.

function make_polynomial_ring(Bs::Vector{Vector{Int}}, B::Vector{Int},
                              F::AbstractAlgebra.Ring)
    
    MC = bases_matrix_coordinates(Bs, B)
    R, x = PolynomialRing(F, :"x"=>MC)
    xdict = Dict{Vector{Int}, MPolyElem}([MC[i] => x[i] for i in 1:length(MC)])
    return R, x, xdict
end



#compute TSC
function TSC_with_basis(M,B,F)
    d = rank(M)
    n = length(matroid_groundset(M))
    Bs = bases(M)
    nBs = nonbases(M)
    R, x, xdict = make_polynomial_ring(Bs,B,F)
    MC = bases_matrix_coordinates(Bs, B)
    
    X = make_coordinate_matrix(d, n, MC,B,R, x,xdict)
    
    Sgens = unique!([det(X[:,b]) for b in Bs])
    Igens = unique!([det(X[:,b]) for b in nBs])
    
    return(Igens,Sgens,R,gens(R),x,xdict)
    
end



# M and N have same number of rows, and M has #B columns, both have entries in ring R
function interlace_columns(M, N, B::Vector{Int},
                           R::MPolyRing, x::Vector{T}) where T <: MPolyElem 
    
    M_nrows, M_ncols = size(M)
    N_nrows, N_ncols = size(N)
    n = M_ncols + N_ncols

    S = MatrixSpace(R, M_nrows, n)
    Bc = [i for i in 1:n if !(i in B)]
    
    X = S()
    X[:, B] = M
    X[:, Bc] = N
    
    return X 
end

# This makes the matrix X from which we compute the coordinate ring of the matroid
# stratum. It has the identity matrix at columns indexed by B, 0's at locations
# determined by the nonbases of X. 
function make_coordinate_matrix(d::Int, n::Int, MC::Vector{Vector{Int}},
                                B::Vector{Int},
                                R::MPolyRing, x::Vector{T},
                                xdict::Dict{Vector{Int}, MPolyElem}) where T <: MPolyElem
    
    Id = identity_matrix(R,d)
    Xpre = make_coordinate_matrix_no_identity(d, n, MC, R, x, xdict)
    return interlace_columns(Id, Xpre, B, R, x)
end


#Coord matrix wrt basis and field
function TSC_matrix(M,B,F)
    
    d = rank(M)
    n = length(matroid_groundset(M))
    Bs = bases(M)
    nBs = nonbases(M)
    R, x, xdict = make_polynomial_ring(Bs,B,F)
    MC = bases_matrix_coordinates(Bs, B)
    
    X = make_coordinate_matrix(d, n, MC,B,R, x,xdict)
    
    return X,R
    
end


#reduced expression for TSC data
function reduce_TSC(M,F)
    
    B = bases(M)[1]
    
    T = TSC(M,B,F)
    
    I = reduce_ideal_full(T[1], T[2], T[3], T[4], false)
    
    return (I[1],I[2])
    
end


#function that finds basis with simplest vanishing ideal, uses that to compute TSC.

function matroid_to_reduced_TSC_min_basis(Q,F)
    
    Bs = bases(Q)
    
    A = argmin(c -> count_nonbases_disjoint_to_chart(Q, c) , Bs)
    
    RQ = TSC_with_basis(Q,A,F)
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
    
    Igens = gens(stepwise_saturation(ideal(Igens), Sgens))
    
    any([is_unit(g) for g in Igens]) && return ([R(1)], Sgens, A)
        
    reducedData = reduce_ideal_full(Igens, Sgens, R, gens(R), false)

    reducedData isa String && return reducedData

    newI = reducedData[1]
    newS = reducedData[2]
    if length(newI) == 0
        newI = [R(0)]
    end
    
    return (newI, newS, A)
    
end