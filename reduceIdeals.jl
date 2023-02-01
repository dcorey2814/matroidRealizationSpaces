
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

   #println("try to isolate = ", v)
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
    
    #if length(with_x_deg_1) == 0
    #    w = [(gen, x, vars(gen), degree(gen,x)) for gen in Igens] 
       #println("w = ", w)
    #end
    
      
   #println("with_x_deg_1 = ", with_x_deg_1)
    
    (length(with_x_deg_1) != 0) || return "can't isolate" 

    for gen in with_x_deg_1
        t = coefficient_v(x,gen,R)        
       #println("t test = ", t)
        if t in S
        	#println("t inverted = ", t)
            return nn_isolate_v(x, gen, R, S)
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
    
    if R(0) in Sgens
        
        return "Not realizable 0 in Semigroup"
        
    else
       #println("Igens1 = ", Igens)
        Ivars = ideal_vars(Igens); 
        
        #println("Ivars")
        #println("Ivars = ", Ivars)
        
        S =  MPolyPowersOfElement(R , Sgens); 
        for x in Ivars 
            tx = find_solution_x(x, Igens, R, S)
           #println("x, find solution = ", x, ",  ", tx)
            if tx isa String
                continue
            else 
            	 #println("x = ", x)
                Igens_new = n_new_Igens(x, tx, Igens, R, varlist); 
                
               #println("Igens_new = ", Igens_new)
                
                Sgens_new = n_new_Sgens(x, tx, Sgens, R, varlist)
    
                #(R(0) in Sgens) && error("Nonrealizable") 
                
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


function clean(f,R)
    
    if !(f == 0)
    
        if length(f) == 1
            
            return R(coeff(f,1))
            
        else
            
            fFactors = factor(f)
    
            FactorsDict = Dict(fFactors)
    
            cleanf = prod([k^(FactorsDict[k]) for k in keys(FactorsDict) if (length(k)>1||is_constant(k))])
        
            return cleanf
        
        end
        
    else
        
        return f
        
    end
    
end

function matroid_to_reduced_expression(Q, F,k)
    
    charts = [c for c in circuits(Q) if length(c) == rank(Q)+1]
    A = charts[1]
    RQ = matroid_realization_space(Q, A, F)
    R = parent(RQ[1][1])
    Sgens = [s for s in RQ[2] if length(s) <= k]#new 13.1.2023
    I = reduce_ideal_full(RQ[1], Sgens, R, gens(R), false)
    
    
 
    if I isa String
        
        return I
        
     #I[1] = ideal generators, I[2] = subgroup generators   
    else
        
       Iclean = unique!([clean(f,R) for f in I[1]])
    
            
       return (Iclean, I[2])
        
       # return(I[1],I[2])
        
    end
end





