
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


#isolate variable v
function nn_isolate_v(v,f,R,S)

    (degree(f,v) == 1) || return "can't isolate, degree of variable must be 1"
    no_v = [term(f,i) for i in 1:length(f) if !(v in vars(monomial(f,i)))]
    (length(no_v) != 0) || return "can't isolate, variable is a factor"
    den = coefficient_v(v,f,R)
    (den in S) || return "can't isolate, coefficient of variable is not a unit"
    h = R(-1)*sum(no_v)
    return h//den

end

#function that finds generators of ideal that contain x, checks if coefficient is in semigroup, then solves
function find_solution_x(x,Igens,R,S)
    
    with_x_deg_1 = [gen for gen in Igens if (x in vars(gen) && degree(gen,x)==1)]    
    (length(with_x_deg_1) != 0) || return "can't isolate" 

    for gen in with_x_deg_1
        t = coefficient_v(x,gen,R)        
        if t in S
            return nn_isolate_v(x,gen,R,S)
        else 
            return "can't solve for variable"
        end
    end
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
function n_new_Igens(x, tx, Igens, R, varlist)    
    return unique!([sub_v(x, tx, gen, R, varlist) for gen in Igens])
end

#function to produce new subgroup generators
function n_new_Sgens(x, tx, Sgens, R, varlist)
    return unique!([sub_v(x, tx, gen, R, varlist) for gen in Sgens])
end

function reduce_ideal_one_step(Igens, Sgens, R, varlist, fullyReduced)
    Ivars = ideal_vars(Igens); 
    S =  MPolyPowersOfElement(R , Sgens); 
    for x in Ivars 
        tx = find_solution_x(x, Igens, R, S)
        if tx isa String
            continue
        else 
            Igens = n_new_Igens(x,tx,Igens, R, varlist); 
            Sgens = n_new_Sgens(x,tx,Sgens,  R, varlist)
#            (R(0) in Sgen) && error("Nonrealizable") 
            return (Igens, Sgens, R, varlist, fullyReduced)
        end
    end
    return (Igens, Sgens, R, varlist, true)
end


function reduce_ideal_full(Igens, Sgens, R, varlist, fullyReduced = false)
    if !fullyReduced
        (Igens, Sgens, R, varlist, fullyReduced) = reduce_ideal_one_step(Igens, Sgens, R, varlist, fullyReduced)
        return reduce_ideal_full(Igens, Sgens, R, varlist, fullyReduced)
    else
        return (Igens, Sgens, R, varlist, fullyReduced)
    end
end




