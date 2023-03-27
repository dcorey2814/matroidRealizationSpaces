##############
## File handling  ##
##############


function string2Int64Vector(s)
    return map(i -> parse(Int64, i), split(s))
end

function file2SetVectors(fileName)
    return map(s -> string2Int64Vector(s), readlines(fileName))
end

function vec2String(v)
    vecString = map(i -> string(i), v)
    return join(vecString, " ")
end
