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

function is_representable(M::Matroid; char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing)::Bool
    RS = new_realization_space(M, char=char, q=q)
    return RS.representable
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

function projective_identity(d::Int)
    if d == 1
        return ones(Int, 1, 1)
    end

    X = zeros(Int, d, d+1)
    for i in 1:d
        X[i,i] = 1
        X[i,d+1] = 1
    end
    return X
end

function interlace_columns(M::MatrixElem{<:RingElem}, N::MatrixElem{<:RingElem}, 
                           B::Vector{Int}, R::AbstractAlgebra.Ring)   
    
    M_nrows, M_ncols = size(M)
    N_nrows, N_ncols = size(N)
    n = M_ncols + N_ncols

    Bc = [i for i in 1:n if !(i in B)]
    
    X = zero_matrix(R, M_nrows, n)
    X[:, B] = M
    X[:, Bc] = N
    
    return X 
end

function realization_matrix(M::Matroid, A::Vector{Int}, F::AbstractAlgebra.Ring)
    
    
    d = rank(M)
    n = length(matroid_groundset(M))
    Bs = bases(M)
    RBC = realization_bases_coordinates(Bs, A)
    D = partial_matrix_max_rows(RBC)
    
    numVars = length(RBC) - (n-d-1)
    R, x = polynomial_ring(F, numVars)
    Id = matrix(R, projective_identity(d))
    
    MR = [x for x in RBC if x[1] != D[x[2]]]
    X = zero_matrix(R, d, n-d-1)
    
    k=1
    for j in 1:n-d-1, i in 1:d
        if [i,j] in MR
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

function new_matroid_realization_space(M::Matroid, A::Vector{Int};
    F::AbstractAlgebra.Ring = ZZ, saturate::Bool=false)::MatroidRealizationSpace

    rk = rank(M)
    n = length(M)
    Bs = bases(M)
    
    polyR, x, mat = realization_matrix(M, A, F)
    
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


#returns the factors of f, but not the exponents
function poly_2_factors(f::RingElem)
    return collect(keys(Dict(factor(f))))
end

# returns the unique factors of the elements of Sgen, again no exponents. 
function gens_2_factors(Sgens::Vector{<:RingElem})
    return unique!(vcat([poly_2_factors(f) for f in Sgens]...))
end


function stepwise_saturation(I::MPolyIdeal, Sgens::Vector{<:RingElem})
    foreach(f -> I = saturation(I,ideal([f])), Sgens)
    return I
end

#####################
# initial reduction #
#####################


function update_hom(v_elim::RingElem, w_repl::RingElem, phi::Oscar.MPolyAnyMap)
    R = codomain(phi)
    phi2_im = replace(gens(R), v_elim => w_repl)
    phi2 = hom(R, R, a->a, phi2_im)
    return hom(R, R, a->a, phi2.(phi.(gens(R))))
end

function small_find_solution(v::RingElem, Igens::Vector{<:RingElem})
    
#    v_deg_1_no_coef = [g for g in Igens if isone(degree(g,v)) && is_unit(coeff(g,[v],[1])) ]
#    println("v_deg_1_no_coef: ", v_deg_1_no_coef)
#    length(v_deg_1_no_coef) != 0 || return "can't isolate"
#    f = first(v_deg_1_no_coef)
#    return -coeff(f, [v], [0]) / coeff(f, [v], [1])
        
    with_v_deg_1 = [g for g in Igens if (total_degree(g) == 1 && is_unit(coeff(g,[v],[1])))] 
    length(with_v_deg_1) == 0 && return "can't isolate"
    
    f = first(with_v_deg_1)
    return -coeff(f, [v], [0])/coeff(f,[v],[1])
end


function small_reduce_one_step(Igens::Vector{<:RingElem}, phi::Oscar.MPolyAnyMap, 
                               elim::Vector{<:RingElem}, fullyReduced::Bool)
    Ivars = ideal_vars(Igens); 
    length(Igens) == 0 && return (Igens, phi, elim, true)
    for v in Ivars 
        w = small_find_solution(v, Igens)
        w isa String && continue
        
        push!(elim, v)
        phi = update_hom(v, w, phi)
        Igens = filter(x->x!=0, phi.(Igens))
        return (Igens, phi, elim, false)
    end
    return (Igens, phi, elim, true)
end


function small_reduce_full_rec(Igens::Vector{<:RingElem}, 
                               phi::Union{<:Oscar.MPolyAnyMap,Nothing} = nothing, 
                               elim::Vector{<:RingElem} = Vector{RingElem}(), 
                               fullyReduced=false)

    if isnothing(phi)
        R = parent(Igens[1])
        phi=hom(R, R, a->a, gens(R))
    end

    output = small_reduce_one_step(Igens, phi, elim, fullyReduced)
    output isa String && return "Not Realizable 0 in Semigroup"
    
    (Igens, phi, elim, fullyReduced) = output    
    !fullyReduced && return small_reduce_full_rec(Igens, phi, elim, fullyReduced)
    
    R = domain(phi)
    x = gens(R)    
    cR = coefficient_ring(R)
    
    if length(elim) == length(x)
        S = cR
    else
        S, y = polynomial_ring(cR, length(x) - length(elim))
    end        
    z=[]
    j=1
    for ele in x
        if ele in elim
            push!(z,0)
        else
            push!(z, y[j])
            j+=1
        end
    end
    phi2 = hom(R, S, a->a, z)
    phif = phi*phi2

    return phif
end


function small_reduce(MRS::MatroidRealizationSpace) 
    
    !MRS.representable && return MRS    
    phi = small_reduce_full_rec(gens(MRS.defining_ideal))
    new_ring = codomain(phi)
    new_ideal = phi(MRS.defining_ideal)
    
    new_ideal_gens = filter!(x->!iszero(x), phi.(gens(MRS.defining_ideal)))
    
    if length(new_ideal_gens) == 0
        new_ideal = ideal(new_ring, [new_ring(0)])
    else
        new_ideal = ideal(new_ideal_gens)
    end

    new_mat = phi.(collect(MRS.representation_matrix))
    new_mat = matrix(new_ring, new_mat)
    new_ineq = unique!(filter(x->!is_unit(x), phi.(MRS.inequations)))
    new_ineq = gens_2_factors(new_ineq)
    
    return MatroidRealizationSpace(new_ideal, new_ineq, new_ring, new_mat, MRS.representable, MRS.F, MRS.char, MRS.q)
    
end




#####################
# full reduction    #
#####################

# computes the coefficient of v in monomial m. 
function coefficient_monomial(v::RingElem, m::RingElem)
    isone(degree(m,v)) || return "The variable is not a degree 1 factor."
    mf = factor(m)
    mfdict = Dict(mf)
    u = unit(mf)
    not_v = [k^(mfdict[k]) for k in keys(mfdict) if k != v ]
    length(not_v) == 0 ? u : u*prod(not_v)
end

# computes the coefficient of v in f. 
function coefficient_v(v::RingElem, f::RingElem)
    isone(degree(f,v)) || return "degree of variable must be 1"
    withv = [term(f,i) for i in 1:length(f) if v in vars(monomial(f,i))]
    return sum([coefficient_monomial(v,m) for m in withv])
end


function find_solution_v(v::RingElem, Igens::Vector{<:RingElem}, 
                         Sgens::Vector{<:RingElem}, R::MPolyRing) 

    
    with_v_deg_1 = [g for g in Igens if isone(degree(g,v))] 
    length(with_v_deg_1) != 0 || return "can't isolate"

    for f in with_v_deg_1

        den = coefficient_v(v, f)
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
