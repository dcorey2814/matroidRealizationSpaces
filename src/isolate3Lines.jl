function lines(Q::Matroid)
    return [l for l in hyperplanes(Q) if length(l) > rank(Q)-1 ]
end

function count_lines_thru_one_point(Ls::Vector{Vector{Int64}}, i::Int64)
    return length([l for l in Ls if i in l])
end

function count_3_lines_thru_all_points(Q::Matroid)
    Ls = lines(Q)    
    return [count_lines_thru_one_point(Ls,i) for i in matroid_groundset(Q)]
end


function to_revlex(Q::Matroid, nCd::Vector{Vector{Int64}})
    B = bases(Q)
    l = []
    for b in nCd
        if b in B
            push!(l,"*")
        else
            push!(l,"0")
        end
    end
    return join(l)
end
