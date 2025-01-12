using SparseArrays, LinearAlgebra
# Compute the boundary of a cubical simplex sx and represent it as a list of maximal simplices, 
# i.e., just determine the codimenstion 1 faces of sx. 
function simplex_boundary(sx)
    borders = []
    for dim in findall([dim % 2 == 1 for dim in sx])
        borders = [borders; [Base.setindex(sx, sx[dim]-1, dim), Base.setindex(sx, sx[dim]+1, dim)]]
    end
    return borders
end


# Create a dictionary of all simplices grouped by dimension for a simplicial complex
# represented by its maximal simplices. I dont think weborders = [borders; [Base.setindex(sx, sx[dim]-1, dim), Base.setindex(sx, sx[dim]+1, dim)]] need this
# function dim_dict_scx(max_scx)
    
# end

function boundary_matrix(dimension_dict, n, p = 0)
    # check if p is 0 or prime
    if p != 0 && !isprime(p)
        throw(DomainError(p, "p should be 0 or a prime."))
    end
    # cases n = 0 and n = dim+1 should be treated separately
    if n == 0
        Dn = zeros(p==0 ? Float64 : Mod{p}, 1, length(dimension_dict[n]))
    elseif n == maximum(keys(dimension_dict))+1
        Dn = zeros(p==0 ? Float64 : Mod{p}, length(dimension_dict[n-1]), 1)
    else
        Dn = zeros(p==0 ? Float64 : Mod{p}, length(dimension_dict[n-1]), length(dimension_dict[n]))
        search = Dict(k => Dict(value => i for (i, value) in enumerate(v)) for (k,v) in dimension_dict)
        # otherwise, just replace the zeroes with 1's or -1's at appropriate places
        for j in eachindex(dimension_dict[n])
            boundary = simplex_boundary(dimension_dict[n][j])
            for i in eachindex(boundary)
                # print((-1)^i)
                # println(boundary[i])
                Dn[search[n-1][boundary[i]], j] = (-1)^(i + floor(i/3))
            end
        end
    end
    return Dn
end

# determine the mod p Betti numbers of a simplicial complex given via its maximal simplices
function betti_numbers(scx, p = 0)
    scx = Dict((x-1)=> scx[x] for x in eachindex(scx))
    # display(scx)
    # determine the dimension of max_scx
    d = maximum(keys(scx))
    # initialize the vector of Betti numbers
    betti = []
    # the current (first) boundary matrix
    Dk1 = boundary_matrix(scx, 0, p)
    rank1 = rank(Dk1)
    # determine the Betti numbers from two consecutive boundary matrices
    for k in 1:(d+1)
        # the next boundary matrix
        Dk = boundary_matrix(scx, k, p)
        # the commented two lines below work only over R (for p = 0)
        # and need 'using LinearAlgebra' for rank computation
        # println("$n_rows, $(Dk.size), $(Dk1.size)")
        # println("$n_rows, $(rank(Dk)), $(rank(Dk1))")
        # perform simultaneous reduction on D_{k-1} and D_k 
        # and store the corresponding Betti number
        # display("Ranks")
        # display(size(Dk, 1))
        # display(rank(Dk1))
        # display(rank(Dk1, atol=1e-10))
        # display(rank(Dk))
        # display(rank(Dk, atol=1e-10))
        # display((-Dk1) * Dk)
        # display(Dk1 * Dk)
        # push!(betti, simultaneous_reduce(Dk1, copy(Dk)))
        # display(nullspace(Dk))
        # display(rowspace(Dk1))
        # display(betti[end])
        n_rows, n_cols = size(Dk)
        rank2 = rank(Dk)
        push!(betti, n_rows - rank2 - rank1)
        # use Dk as the current (first) boundary matrix
        rank1 = rank2
        Dk1 = Dk
    end
    return betti
end

# This function performs a simultaneous reduction on matrices A=D_n and B=D_{n+1}
# to determine the number of generators of H_n. It implicitly assumes that AB=0.
function simultaneous_reduce(A, B)
    # throw a 'wrong size' error
    if size(A)[2] != size(B)[1]
        throw(DimensionMismatch("cols(A) != rows(B)"))
    end
    # display(A)
    # display(B)
    # store the number of rows and columns of A 
    n_rows, n_cols = size(A)
    # perform a column reduction on A and do the inverse operations on rows of B
    i, j = 1, 1
    last_nonzero_col = 0
    while i <= n_rows && j <= n_cols
        # check if A[i, j] can be used as a column pivot ...
        if A[i, j] == 0
            # if not, find a column pivot in the current row
            nonzero_col = j+1
            while nonzero_col <= n_cols && A[i, nonzero_col] == 0
                nonzero_col += 1;
            end
            # if there are only zeros in the current row to the right
            if nonzero_col == n_cols+1
                # go to a new row
                i += 1
                continue
            end
            # otherwise swap the columns of A ...
            A[:, [j nonzero_col]] = A[:, [nonzero_col j]]
            # ... and corresponding rows of B
            B[[j nonzero_col], :] = B[[nonzero_col j], :]
            # code here ...
        end
        # store last nonzero column (for additional row reduction on B later)
        last_nonzero_col = j
        # ... and set A[i, j] as the new column pivot entry 
        pivot = A[i, j]
        # set that column pivot to 1 (in A)
        A[:, j] = A[:, j]/pivot
        # and do the inverse operation on B
        B[j, :] = B[j, :]*pivot
        # annihilate the nonzero elements to the right of the pivot in A
        # and do the opposite transformations on B
        for col in (j+1):n_cols
            scale = A[i, col]
            A[:, col] -= scale*A[:, j]
            B[j, :] += scale*B[col, :]
        end
        # increase row and column indices and search for the next pivot
        i += 1
        j += 1
    end
    # display(A)
    # display(B)
    # now perform a row reduction on nonzero rows of B 
    # (Performing inverse operations on columns of A would only be required to keep track of the generators...)
    n_rows, n_cols = size(B)
    i, j = last_nonzero_col+1, 1
    last_nonzero_row = last_nonzero_col
    while i <= n_rows && j <= n_cols
        # check if B[i, j] can be used as a column pivot ...
        if B[i, j] == 0
            # if not, find a row pivot in the current column
            nonzero_row = i+1
            while nonzero_row <= n_rows && B[nonzero_row, j] == 0
                nonzero_row += 1
            end
            # if there are only zeros in the current column downwards
            if nonzero_row == n_rows+1
                # go to a new column
                j += 1
                continue
            end
            # otherwise swap the rows of B
            B[[i nonzero_row], :] = B[[nonzero_row i], :]
            # (The corresponding columns of the column-reduced A are 0, no need to swap those...)
        end
        # store last nonzero row (to return the 'Betti number')
        last_nonzero_row = i
        # Set B[i, j] as the new pivot entry ...
        pivot = B[i, j]
        # ... and set that pivot to 1.
        B[i, :] = B[i, :]/pivot
        # annihilate the nonzero elements below the pivot in B
        # (No need for inverse operations on the columns of the reduced A...)
        for row in (i+1):n_rows
            scale = B[row, j]
            B[row, :] -= scale*B[i, :]
        end
        # increase row and column indices and search for the next pivot
        i += 1
        j += 1
    end
    # display(A)
    # display(B)
    # return the corresponding Betti number
    # println("$n_rows $last_nonzero_row")
    return n_rows - last_nonzero_row
end