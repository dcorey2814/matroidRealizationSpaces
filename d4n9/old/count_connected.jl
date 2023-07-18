d4n9 = vec(readlines("d4n9.dat"))

for a in d4n9 

    if is_connected(matroid_from_revlex_basis_encoding(a,4,9))

                    open("connected_4_9.dat", "a") do file
                         write(file,a,"\n")
                    end 
    end
end

print("done")

