function maximal_circuits(Q)
    r = rank(Q)
    return [c for c in circuits(Q) if length(c) == r+1]
end

function count_nonbases_chart_int2(Q, A)
    NBs = nonbases(Q)
    r = rank(Q)
    return length([nb for nb in NBs if length(intersect(A,nb)) == r-1])
end

function matroid_with_chart_to_reduced_expression(Q, Ac, Fc)     
    MRS = new_matroid_realization_space(Q, Ac; F = Fc)
    return reduce_ideal_full(MRS)
end












