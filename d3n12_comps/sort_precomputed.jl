currentDir = pwd()

d3n12 = vec(readlines(joinpath(currentDir, "d3n12/lines_3.dat")));

not_principal = vec(readlines(joinpath(currentDir, "d3n12/precomputed/not_principal_multivariate_pc.dat")));
not_realizable = vec(readlines(joinpath(currentDir, "d3n12/precomputed/not_realizable_pc.dat")));
princ_multi = vec(readlines(joinpath(currentDir, "d3n12/precomputed/principal_multivariate_pc.dat")));
princ_univ = vec(readlines(joinpath(currentDir, "d3n12/precomputed/principal_univariate_pc.dat")));
zero_ideal = vec(readlines(joinpath(currentDir, "d3n12/precomputed/zero_ideal_pc.dat")));

pc = sort!(vcat(not_principal, not_realizable, princ_multi, princ_univ, zero_ideal));

to_sort = setdiff(d3n12, pc)

io = open(joinpath(currentDir, "d3n12/to_sort.dat"), "w")

for a in to_sort
    write(io, a*"\n")
end

close(io)
