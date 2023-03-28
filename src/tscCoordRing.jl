##########################################
##   Creating the coordinate ring of thin Schubert cells     ###
##########################################

function makePolyRing(d,n,F)
    R,x = PolynomialRing(F, :"x"=>(1:d, 1:n-d))
    W = ones(Int64, d*(n-d))
    R,x = grade(R,W)
    x = reshape(x, (d,n-d))
    return R,x
end


function makeCoordinateMatrix(d,n,B,R,x)
    Id = identity_matrix(R,d)
    S_id = MatrixSpace(R,d,n)
    S = MatrixSpace(R,d,n-d)
    X = S(x)
    
    X_id = S_id()
    
    idCount = 1
    varCount = 1
    
    for i in 1:n
        if i in B
            X_id[:,i] = Id[:,idCount]
            idCount += 1
        end
        if i ∉ B
            X_id[:,i] = X[:,varCount]
            varCount +=1
        end
    end
    
    return  X_id
end

function reducedCoordinateMatrix(M, d, n, B, R, x)
    
    X = makeCoordinateMatrix(d,n,B,R,x)
    nonBasesM = nonbases(M)
    coordNonBases = [ nb for nb in nonBasesM if length(symdiff(B,nb)) == 2]
    
    if length(coordNonBases) == 0
        return X
    end
    
    for nb in coordNonBases
        row_nb = setdiff(B,nb)[1]
        row_nb = length([a for a in B if a < row_nb]) + 1
        
        col_nb = setdiff(nb,B)[1]
        X[row_nb,col_nb] = 0
    end
    return X
    
end

function TSC(M, F, B, R, x)
    d = rank(M)
    n = length(matroid_groundset(M))
    Bs = bases(M)
    NBs = nonbases(M)
    #R,x = PolynomialRing(F, :"x"=>(1:d, 1:n-d))

    Mx = reducedCoordinateMatrix(M,d,n,B,R,x)
    I = ideal(unique!([det(Mx[:, nb]) for nb in NBs ]))
    toInvert = unique!([det(Mx[:, b])  for b in Bs])
    
    for s in toInvert
        I = saturation(I,ideal([s]))
    end
    
    return I
end

function projectedNonBases(M, F, B, R, x)
    d = rank(M)
    n = length(matroid_groundset(M))
    NBs = nonbases(M)
    Mx = reducedCoordinateMatrix(M,d,n,B,R,x)
    return unique!([det(Mx[:, nb]) for nb in NBs ]) 
end


function projectedBases(M, F, B, R, x)
    d = rank(M)
    n = length(matroid_groundset(M))
    Bs = bases(M)
    Mx = reducedCoordinateMatrix(M,d,n,B,R,x)
    return unique!([det(Mx[:, nb]) for nb in Bs ])
end

function limitTSC(Ms, F, B, R, x)
    d = rank(Ms[1])
    n = length(matroid_groundset(Ms[1]))

    nonBasesX = []
    basesX = []

    for M in Ms        
        append!(nonBasesX, projectedNonBases(M, F, B, R, x))
        append!(basesX, projectedBases(M, F, B, R, x))
    end
    
    I = ideal(R, unique!(nonBasesX))
    
    
    for s in basesX
        I = saturation(I,ideal([s])); 
    end
    
    return I
end
#new TSC function
function new_TSC(M, F)
    d = rank(M)
    n = length(matroid_groundset(M))
    Bs = bases(M)
    NBs = nonbases(M)
    B = Bs[1]
    R,x = makePolyRing(d,n,F)
   
    Mx = reducedCoordinateMatrix(M,d,n,B,R,x)
    I = ideal(unique!([det(Mx[:, nb]) for nb in NBs ]))
    toInvert = unique!([det(Mx[:, b])  for b in Bs])
    
    return (I,toInvert)
end
#TSC ideal reduction