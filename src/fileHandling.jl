function string2Int64Vector(s::String)
    return map(i -> parse(Int64, i), split(s))
end

function vec2String(v::Vector{Vector{Int64}})
    vecString = map(i -> string(i), v)
    return join(vecString, " ")
end

function file2SetVectors(fileName::String)
    return map(s -> string2Int64Vector(s), readlines(fileName))
end

function to_star0(S::Vector, n::Int64)
    return join(map(x -> x in S ? "*" : "0", 1:n ))
end
