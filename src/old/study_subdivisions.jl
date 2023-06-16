#code to compute and study corank subdivisions of hyperimplex

#map vertex of polytope to indices of nonzero entries
function vertex_to_support(x)
    set = []
    for t in 1:length(x)
        if !(x[t] == 0)
            push!(set,t)
        end
    end
    return set
end

#compute corank vector of a matroid
function corank_vector(M)
    d = rank(M)
    n = length(matroid_groundset(M))
    
    D = pm.polytope.hypersimplex(d,n)
    V = D.VERTICES[:,[i for i in 2:n+1]]
    
    
    x = zeros(length(V[:,1]))

    
    for t in 1:length(V[:,1])
        
        b = vertex_to_support(V[t,:])
        
        x[t] = d - rank(M,b)#corank
        
        
    end
    return x
end

#checks if polytope is generalized permutahedron. input polytope needs to be convex_hull(point configuration)

function is_genperm(P)
    
    E = [vertices(e) for e in faces(P,1)]
    
    yes = []
    
    for edge in E
        
        diff = edge[1]-edge[2]
        
        supp = [i for i in 1:length(diff) if !(diff[i] == 0)]
        
        if !(length(supp) == 2)
            
            return 1 == 2
            
            
        else
            
             if diff[supp[1]] == -diff[supp[2]]
                
                push!(yes,edge)
                
            end
            
        end
        
    end

     return length(E) == length(yes)
    
end
#this is way too slow and clunky, but it works


#input: matroid subdivision of hypersimplex, rank and dimension. returns matroids corresponding to cells.
function subdivision_to_matroids(S,d,n)
    matroids = []
    
    Cells = maximal_cells(S)
    D = pm.polytope.hypersimplex(d,n)
    
    
    V = D.VERTICES[:,[n for n in 2:n+1]]
    
    for c in Cells
        verts = [vertex_to_support(V[t,:]) for t in c]
        M = matroid_from_bases(verts,n)
        push!(matroids,M)
    end
    
    return matroids
    
end