# Functions in this file compute the homology groups of a simplicial complex
# with Z/pZ coefficients for a prime p or R (real) coefficients
# Loading of 'Primes' and 'Mods' packages is assumed via 'using Primes, Mods' in the Julia REPL.

# Compute the boundary of a simplex sx and represent it as a list of maximal simplices, 
# i.e., just determine the codimenstion 1 faces of sx. 
# E.g. for sx = (1, 2, 3) we get [(1, 2), (1, 3), (2, 3)].
function simplex_boundary(sx)
    # determine the codimension 1 simplices in the boundary of the simplex sx
    # 'remove' one entry of a tuple at each specific index by splicing the corresponding subtuples
    cd1faces = ([(sx[1:k-1]..., sx[k+1:end]...) for k in eachindex(sx)])
    return cd1faces
end

# Create a dictionary of all simplices grouped by dimension for a simplicial complex
# represented by its maximal simplices.
function dim_dict_scx(max_scx)
    # determine the dimension+1 of the complex
    n = maximum([length(sx) for sx in max_scx])
    # set up the dictionary (empty vectors of k-tuples in each dimension)
    dimension_dict = Dict((k-1) => Vector{NTuple{k, Int}}() for k in 1:n)
    # add all simplices in max_scx to the dictionary
    # (Sorting of tuple entries is needed to ease recognition of generators later.)
    for sx in max_scx
        push!(dimension_dict[length(sx)-1], Tuple(sort(collect(sx))))
    end
    # fill the dictionary with all simplices from top down
    for k in (n-1):-1:1
        for sx in dimension_dict[k]
            union!(dimension_dict[k-1], simplex_boundary(sx))
        end
    end
    return dimension_dict
end

# This function returns the n-th boundary matrix Dn with coefficients Z/pZ or R of a simplicial complex given as a dimension_dict.
# (p should be a prime. For convenience, p = 0 means R.)
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
        # otherwise, just replace the zeroes with 1's or -1's at appropriate places
        for j in eachindex(dimension_dict[n])
            boundary = simplex_boundary(dimension_dict[n][j])
            for i in eachindex(boundary)
                Dn[findfirst(==(boundary[i]), dimension_dict[n-1]), j] = (-1)^i
            end
        end
    end
    return Dn
end

# determine the mod p Betti numbers of a simplicial complex given via its maximal simplices
function betti_numbers(max_scx, p = 0)
    # build the 'dimension dictionary' of a simplicial complex given via its maximal simplices
    dimension_dict = dim_dict_scx(max_scx)
    # determine the dimension of max_scx
    d = maximum(keys(dimension_dict))
    # initialize the vector of Betti numbers
    betti = []
    # the current (first) boundary matrix
    Dk1 = boundary_matrix(dimension_dict, 0, p)
    # determine the Betti numbers from two consecutive boundary matrices
    for k in 1:(d+1)
        # the next boundary matrix
        Dk = boundary_matrix(dimension_dict, k, p)
        # the commented two lines below work only over R (for p = 0)
        # and need 'using LinearAlgebra' for rank computation
        #n_rows, n_cols = size(Dk)
        #push!(betti, n_rows - rank(Dk) - rank(Dk1))
        # perform simultaneous reduction on D_{k-1} and D_k 
        # and store the corresponding Betti number
        push!(betti, simultaneous_reduce(Dk1, Dk))
        # use Dk as the current (first) boundary matrix
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
    # return the corresponding Betti number
    return n_rows - last_nonzero_row
end