#code that computes the coordinate ring of a TSC with respect to some matroid, basis, and field. We can also give the associated matrix with respect to the same basis and field. Also, we can give the reduced presentation of the coordinate ring given some basis and field.

#Note that this code needs to interact with the code from matroid_realization.jl

#Throughout M = matroid, B = basis, F= field.

################################################
#compute TSC
function TSC(M,B,F)
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


#Coord matrix wrt basis and field
function TSC_matrix(M,B,F)
    
    d = rank(M)
    n = length(matroid_groundset(M))
    Bs = bases(M)
    nBs = nonbases(M)
    R, x, xdict = make_polynomial_ring(Bs,B,F)
    MC = bases_matrix_coordinates(Bs, B)
    
    X = make_coordinate_matrix(d, n, MC,B,R, x,xdict)
    
    return X
    
end


#reduced expression for TSC data
function reduce_TSC(M,F)
    
    B = bases(M)[1]
    
    T = TSC(M,B,F)
    
    I = reduce_ideal_full(T[1], T[2], T[3], T[4], false)
    
    return (I[1],I[2])
    
end