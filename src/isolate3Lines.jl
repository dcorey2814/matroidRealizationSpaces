function lines(M::Matroid)
    return [l for l in hyperplanes(M) if length(l) > rank(M)-1 ]
end

function count_lines_thru_one_point(Ls,i)
    return length([l for l in Ls if i in l])
end

function count_3_lines_thru_all_points(M)
    Ls = lines(M)    
    return [count_lines_thru_one_point(Ls,i) for i in matroid_groundset(M)]
end



function to_revlex(M, nCd)
    B = bases(M)
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
