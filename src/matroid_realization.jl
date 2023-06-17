struct MatroidRealizationSpace
    defining_ideal::Union{Ideal, NumFieldOrdIdl, Nothing}
    inequations::Union{Vector{Oscar.RingElem},Nothing}
    ambient_ring::Union{Oscar.MPolyRing, Oscar.Ring, Nothing}
    representation_matrix::Union{Oscar.MatElem,Nothing}
    representable::Bool
end

function Base.show(io::IO, RS::MatroidRealizationSpace)
    if !RS.representable
        print(io, "The matroid is not representable.")
    else
        println(io, "The representations of the matroid are parametrized by the matrix")
        # println isn't ideal as it prints the matrix as one big line
        display(RS.representation_matrix)
        println(io, "in the ", RS.ambient_ring)
        println(io, "within the vanishing set of the ideal\n", RS.defining_ideal)
        println(io, "avoiding the zero loci of the polynomials\n", RS.inequations)
    end
end

#returns the factors of f, but not the exponents
function poly_2_factors(f::RingElem)
    return collect(keys(Dict(factor(f))))
end

# returns the unique factors of the elements of Sgen, again no exponents. 
function gens_2_factors(Sgens::Vector{<:RingElem})
    return unique!(vcat([poly_2_factors(f) for f in Sgens]...))
end


function stepwise_saturation(I::Ideal, Sgens::Vector{<:RingElem})
    foreach(f -> I = saturation(I,ideal([f])), Sgens)
    return I
end



function projective_identity(d::Int)
    if d == 1
        return ones(Int, 1, 1)
    end

    X = zeros(Int, d, d+1)
    for i in 1:d
        X[i, i] = 1
        X[i, d+1] = 1
    end
    return X
end

function interlace_columns(M::MatrixElem{<:RingElem}, N::MatrixElem{<:RingElem}, 
                           A::Vector{Int}, R::AbstractAlgebra.Ring)   
    
    M_nrows, M_ncols = size(M)
    N_nrows, N_ncols = size(N)
    n = M_ncols + N_ncols

    Ac = [i for i in 1:n if !(i in A)]
    
    X = zero_matrix(R, M_nrows, n)
    X[:, A] = M
    X[:, Ac] = N
    
    return X 
end

function realization_bases_coordinates(Bs::Vector{Vector{Int}}, A::Vector{Int})

    d = length(Bs[1])
    B = A[1:d]
    c1 = A[d+1]

    coord_bases = [b for b in Bs if length(symdiff(B,b)) == 2]    
    new_coords = Vector{Vector{Int}}()

    for b in coord_bases
        if is_subset(b, A)
            continue
        end
        
        row_b = setdiff(B,b)[1]
        row_b = count(a->(a<row_b), B) + 1
        
        col_b = setdiff(b,B)[1]
        col_b = col_b - length([a for a in A if a <= col_b])

        push!(new_coords, [row_b, col_b])
    end
    
    return sort!(new_coords, by = x -> (x[2], x[1]))
end


function partial_matrix_max_rows(Vs::Vector{Vector{Int}})
    nr = maximum([x[1] for x in Vs])
    cols = unique!([x[2] for x in Vs])
    first_nonzero_cols = Dict{Int, Int}(c => maximum(i for i in 1:nr if [i,c] in Vs) for c in cols)
    return first_nonzero_cols 
end


function realization_matrix(Q::Matroid, A::Vector{Int}, F::AbstractAlgebra.Ring)
    d = rank(Q)
    n = length(matroid_groundset(Q))
    Bs = bases(Q)
    RBC = realization_bases_coordinates(Bs, A)
    D = partial_matrix_max_rows(RBC)
    
    numVars = length(RBC) - (n-d-1)
    R, x = polynomial_ring(F, numVars)
    Id = matrix(R, projective_identity(d))
    
    QR = [x for x in RBC if x[1] != D[x[2]]]
    X = zero_matrix(R, d, n-d-1)
    
    k=1
    for j in 1:n-d-1, i in 1:d
        if [i,j] in QR
            X[i,j] = x[k]; k+=1
        elseif(j in keys(D) && i == D[j])
            X[i,j] = R(1)            
        else
            X[i,j] = R(0)
        end
    end
    
    mat = interlace_columns(Id, X, A, R)
    return (R, x, mat)
end

function new_matroid_realization_space(Q::Matroid, A::Vector{Int};
    F::AbstractAlgebra.Ring = ZZ, saturate::Bool=false)::MatroidRealizationSpace

    rk = rank(Q)
    n = length(Q)
    Bs = bases(Q)
    
    polyR, x, mat = realization_matrix(Q, A, F)
    
    eqs = Vector{RingElem}()
    ineqs = Vector{RingElem}()
    
    
    for col in subsets(Vector(1:n),rk)    
        col_det = det(mat[:,col])
        
        if total_degree( col_det ) <= 0 

            if col_det != 0 && col in Bs 
                isunit(col_det) && continue
            elseif col_det != 0 # and col is not a basis
                error("determinant nonzero but set not a basis")
            elseif col in Bs 
                error( "determinant zero but set is a basis" )
            else
                continue
            end
        end
        
        if  col in Bs
            push!(ineqs, col_det)
        else
            push!(eqs, col_det)
        end
    end
    
    def_ideal = ideal(polyR, eqs)
    def_ideal = ideal(groebner_basis(def_ideal))
    ineqs = gens_2_factors(ineqs)
    
    if saturate #|| polyR.nvars < 10
        def_ideal = stepwise_saturation(def_ideal,ineqs)
        def_ideal = ideal(groebner_basis(def_ideal))
    end
    
    representable = !(isone(def_ideal))

    !representable && return MatroidRealizationSpace(def_ideal, nothing, nothing, nothing, false)
    
    return MatroidRealizationSpace(def_ideal, ineqs, polyR, mat, representable)
    
end


function rank_plus1_circuits(Q::Matroid)
    rk = rank(Q) 
    return [c for c in circuits(Q) if length(c) == rk+1]
end

function count_nonbases_meeting_circuit_2elem(Q::Matroid, A::Vector{Int64})
    NBs = nonbases(Q)
    return length([nb for nb in NBs if length(intersect(A,nb)) == 2])
end

function optimal_circuits(Q::Matroid)
    As = rank_plus1_circuits(Q)
    mx = maximum([count_nonbases_meeting_circuit_2elem(Q, A) for A in As])
    return [A for A in As if count_nonbases_meeting_circuit_2elem(Q, A) == mx] 
end

#####################
# full reduction    #
#####################


function find_solution_v(v::RingElem, Igens::Vector{<:RingElem}, 
                         Sgens::Vector{<:RingElem}, R::MPolyRing) 

    
    with_v_deg_1 = [g for g in Igens if isone(degree(g,v))] 
    length(with_v_deg_1) != 0 || return "can't isolate"

    for f in with_v_deg_1

        den = coeff(f, [v], [1])
        fac_den = poly_2_factors(den)
        !issubset(fac_den, Sgens) && continue
	
	
        no_v = [term(f,i) for i in 1:length(f) if !(v in vars(monomial(f,i)))]
        iszero(length(no_v)) && continue
              
        h = R(-1)*sum(no_v)
        return h//den
        
    end
    return "can't solve for v"
end


# v is replaced by t in f
function sub_map(v::RingElem, t::RingElem, R::MPolyRing, xs::Vector{<:RingElem}) 
    xs_v = map(x -> x==v ? t : x, xs )    
    return hom(R,FractionField(R), a->a, xs_v)
end



# replace v by t in f, only return the numerator.
function sub_v(v::RingElem, t::RingElem, f::RingElem, R::AbstractAlgebra.Ring, xs::Vector{<:RingElem}) 
    m = sub_map(v,t,R,xs) 
    new_f = numerator(m(f))
    return new_f 
end


# removes factors that are in the semigroup generated by Sgens
function clean(f::RingElem, R::MPolyRing, Sgens::Vector{<:RingElem})   
    
    fFactors = factor(f)
    FactorsDict = Dict(fFactors) 
    cleanf_arr = [k^(FactorsDict[k]) for k in keys(FactorsDict) if !(k in Sgens) || is_unit(k)]
    
    length(cleanf_arr) > 0 ? prod(cleanf_arr) : unit(fFactors)
    
end

#variables in ideal
function ideal_vars(Igens::Vector{<:RingElem}) 
    return unique!(vcat([vars(gen) for gen in Igens]...))
end

function n_new_Sgens(x::RingElem, t::RingElem, Sgens::Vector{<:RingElem}, 
                     R::AbstractAlgebra.Ring, xs::Vector{<:RingElem}) 
    preSgens = unique!([sub_v(x, t, f, R, xs) for f in Sgens])
    return gens_2_factors(preSgens)
end

function n_new_Igens(x::RingElem, t::RingElem, Igens::Vector{<:RingElem}, 
                     Sgens::Vector{<:RingElem}, R::AbstractAlgebra.Ring, xs::Vector{<:RingElem}) 

    preIgens = unique!([clean(sub_v(x, t, f, R, xs), R, Sgens) for f in Igens])
    return filter(x-> x!= R(0), preIgens)
end



function matrix_clear_den_in_col(X::Oscar.MatElem, c::Int)

    Xc = [denominator(f) for f in X[:, c]]
    t = lcm(Xc)
    return multiply_column!(X, t, c)

end



function matrix_clear_den(X::Oscar.MatElem)
    rs, cs = size(X)
    for c in 1:cs
        X = matrix_clear_den_in_col(X, c)
    end
    return X
end


function reduce_ideal_one_step(MRS::MatroidRealizationSpace, 
                               elim::Vector{<:RingElem}, 
                               fullyReduced::Bool)

    Igens = gens(MRS.defining_ideal)
    Sgens = MRS.inequations
    R = MRS.ambient_ring
    FR = FractionField(R)
    xs = gens(R)
    X = MRS.representation_matrix
    nr, nc = size(X)
    
    Ivars = ideal_vars(Igens);

    for x in Ivars 
        t = find_solution_v(x, Igens, Sgens, R)
        t isa String && continue
        
        phi = sub_map(x, t, R, xs)
        
	Sgens_new = n_new_Sgens(x, t, Sgens, R, xs);
        Igens_new = n_new_Igens(x, t, Igens, Sgens_new, R, xs);
        push!(elim, x)
        
        phiX = matrix(FR, [phi(X[i,j]) for i in 1:nr, j in 1:nc  ] )
        nX_FR = matrix_clear_den(phiX)
        nX = matrix(R, [numerator(nX_FR[i,j])  for i in 1:nr, j in 1:nc ])
        
        GBnew = collect(groebner_basis(ideal(R, Igens_new)))         
        
        MRS_new = MatroidRealizationSpace(ideal(R, GBnew), Sgens_new, R, nX, MRS.representable)

        
        return (MRS_new, elim, fullyReduced)
    end

    return (MRS, elim, true)

end


function reduce_ideal_full(MRS::MatroidRealizationSpace,
                               elim::Vector{RingElem} = Vector{RingElem}(), 
                               fullyReduced::Bool = false) 
        
    output = reduce_ideal_one_step(MRS, elim, fullyReduced)
    output isa String && return "Not Realizable 0 in Semigroup"
    (MRS, elim, fullyReduced) = output    
    
    !fullyReduced && return reduce_ideal_full(MRS, elim, fullyReduced)
    
    
    R = MRS.ambient_ring
    xs = gens(R)
    cR = coefficient_ring(R)
    X = MRS.representation_matrix
    nr, nc = size(X)
    Igens = gens(MRS.defining_ideal)
    Sgens = MRS.inequations

    
    zero_elim = []        
    for i in 1:length(xs)
        if xs[i] in elim
            push!(zero_elim, 0)
        else
            push!(zero_elim, "x"*string(i) ) 
        end
    end
        
    xnew_str = Vector{String}(filter(x -> x!=0,  zero_elim))    
        
    if length(xnew_str) == 0
        phi = hom(R, cR, a->a, [cR(0) for i in 1:length(xs)])
    else
        Rnew, xnew = polynomial_ring(coefficient_ring(R), length(xnew_str)) 
    
        zero_elim_var = []
        j=1
        for i in 1:length(zero_elim)
            if xs[i] in elim
                push!(zero_elim_var, Rnew(0))
            else
                push!(zero_elim_var, xnew[j] ) 
                j+=1
            end
        end
        
        phi = hom(R, Rnew, a->a, zero_elim_var)
    
    end
    
    ambR = codomain(phi)
    Inew = ideal(ambR, phi.(Igens))
    Sgens_new = phi.(Sgens)
    
    normal_Sgens = gens_2_factors([normal_form(g, Inew) for g in Sgens_new])
    unique!(normal_Sgens)


    Xnew = matrix(ambR, [phi(X[i,j]) for i in 1:nr, j in 1:nc])

    MRS_new = MatroidRealizationSpace(Inew, normal_Sgens, ambR, Xnew, MRS.representable)

    return MRS_new
end
