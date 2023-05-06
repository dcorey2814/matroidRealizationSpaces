
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


#remove unit factors
function clean(f, R, Sgens)
    #println(f, " clean enter")
    
    fFactors = factor(f)
    FactorsDict = Dict(fFactors)
    
    cleanf_arr = [k^(FactorsDict[k]) for k in keys(FactorsDict) if !(k in Sgens) || is_constant(k)]
    
    if length(cleanf_arr) > 0
    	cleanf = prod(cleanf_arr)
    else
        cleanf = unit(fFactors)
    end
    #println(cleanf, " clean exit")
    return cleanf 
end
 
#data for realization spaces after reduction


function count_nonbases_disjoint_to_chart(Q, A)
    NBs = nonbases(Q)
    return length([nb for nb in NBs if length(intersect(A,nb)) == 0])
end


function count_nonbases_chart_int2(Q, A)
    NBs = nonbases(Q)
    return length([nb for nb in NBs if length(intersect(A,nb)) == 2])
end



function stepwise_saturation(I, Sgens)
    for f in Sgens
        I = saturation(I, ideal([f]))
    end
    return I
end

function matroid_with_chart_to_reduced_expression(Q, A, F)
     
    RQ = matroid_realization_space(Q, A, F)
    R = parent(RQ[1][1])
    Sgens = RQ[2]    
    Sgens = gens_2_factors(Sgens)
    
    Igens_notsat = gens(ideal(RQ[1]))
    reducedData = reduce_ideal_full(Igens_notsat, Sgens, R, gens(R), false)
    
    if reducedData isa String
        return reducedData
    end
    
    Igens = reducedData[1]
    Sgens = reducedData[2]
    
    
    if length(Igens)==0
        Igens = [R(0)]
    end
    
    Igens = gens(stepwise_saturation(ideal(Igens), Sgens))
    
    if R(1) in Igens
        return ([R(1)], Sgens, A)
    end    
    
    reducedData = reduce_ideal_full(Igens, Sgens, R, gens(R), false)

    if reducedData isa String
        return reducedData
    else    

        newI = reducedData[1]
        newS = reducedData[2]
        if length(newI) == 0
            newI = [R(0)]
        end
        return (newI, newS, A)
    end
end
    
    
#    Igens = gens(stepwise_saturation(ideal(RQ[1]), Sgens))
    
#    if R(1) in Igens
#        return ([R(1)], Sgens, A)
#    end    
    
#    reducedData = reduce_ideal_full(Igens, Sgens, R, gens(R), false)
     
#    if reducedData isa String
#        return reducedData
#    else
       
#       newI = reducedData[1]
#       newS = reducedData[2]
       
#       return (newI, newS, A)
#    end

function maximal_circuits(Q)
    r = rank(Q)
    return [c for c in circuits(Q) if length(c) == r+1]
end


function matroid_to_reduced_expression(Q, F)
    
    charts = [c for c in circuits(Q) if length(c) == rank(Q)+1]    
    overlapCharts = [count_nonbases_chart_int2(Q, c) for c in charts]
    A = argmax(c -> count_nonbases_chart_int2(Q, c) , charts)
    ISA = matroid_with_chart_to_reduced_expression(Q, A, F)
    
    return (ISA[1], ISA[2]) 
     
end















function matroid_to_reduced_expression_min_chart(Q, F)
    
    charts = [c for c in circuits(Q) if length(c) == rank(Q)+1]    
    A = argmin(c -> count_nonbases_disjoint_to_chart(Q, c) , charts)
    ISA = matroid_with_chart_to_reduced_expression(Q, A, F)
    
    return (ISA[1], ISA[2]) 
     
end


function matroid_to_reduced_expression_fewest_var(Q, F)
    
    charts = [c for c in circuits(Q) if length(c) == rank(Q)+1]    
    n = maximum([count_nonbases_chart_int2(Q, c) for c in charts])
    As = [c for c in charts if count_nonbases_chart_int2(Q, c) >= n-1]    
    
    IS = []
    for A in As
        to_add = matroid_with_chart_to_reduced_expression(Q, A, F)
        if length(to_add[1]) == 0
            return to_add
        end
        push!(IS, to_add)
    end
#    IS = [matroid_with_chart_to_reduced_expression(Q, A, F) for A in As]
    
        
    return argmin(a -> (ideal_vars(a[1]), length(a[1]), sum([length(b) for b in a[1]]) ) , IS)
     
end

#reduce TSC ideals
function TSC_to_reduced_expression(M, F)
    
    charts = bases(M)
    #println("1")
    A = charts[1]
    #println("2")
    RQ = new_TSC(M,F)
    #println("3")
    R = parent(RQ[1][1])
    #println("4")
    Sgens = [s for s in RQ[2]]#new 13.1.2023
    #print(length(Sgens),"\n")
    #print(length(gens(RQ[1])))
    #println("5")
    I = reduce_ideal_full(gens(RQ[1]), Sgens, R, gens(R), false)
    #println("6")
    
     varlist = 
    if I isa String
        
        return I
        
     #I[1] = ideal generators, I[2] = subgroup generators   
    else
        
       Iclean = unique!([clean(f,R) for f in I[1]])
    
            
       return (Iclean, I[2])
        
       # return(I[1],I[2])
        
    end
end
