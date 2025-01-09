# Compute the boundary of a cubical simplex sx and represent it as a list of maximal simplices, 
# i.e., just determine the codimenstion 1 faces of sx. 
function same_axis(a, b)
    return sum(map(x -> if x[1] == x[2] 1 else 0 end, zip(a, b)))
end

function simplex_boundary(sx)
    sx = sx |> collect
    if length(sx) == 1
        return Set(sx)
    elseif length(sx) == 2
        return Set(sx)
    elseif length(sx) == 4
        one_corner = sx[1]
        other_corner = filter(x -> same_axis(x, one_corner) == 1, sx |> collect)[1]
        sxes = []
        for corner in [one_corner, other_corner]
            append!(sxes, map(x -> (corner, x),  filter(x -> same_axis(x, corner) == 2, sx)))
        end
        return Set(map(x -> Set(x), sxes))
    elseif length(sx) == 8
        # println(sx)
        one_corner = sx[1]
        other_corner = filter(x -> same_axis(x, one_corner) == 0, sx)[1]
        sxes = []
        for corner in [one_corner, other_corner]
            for i in 1:3
                push!(sxes, filter(x -> x[i] == corner[i], sx))
            end
        end
        # println(Set(map(x -> Set(x), sxes)))
        return Set(map(x -> Set(x), sxes))
    else
        println(sx)
        println("makaj")
        return
    end
end

function to_sets(scx)
    (vertices, edges, squares, cubes) = scx
    vertices = Set(vertices)
    edges = Set(map(x -> Set(x), edges))
    squares = Set(map(x -> Set(x), squares))
    cubes = Set(map(x -> Set(x), cubes))
    return (vertices, edges, squares, cubes)
end

function maximal_sxes(scx::Tuple)
    (vertices, edges, squares, cubes) = scx |> to_sets
    cube_bound = reduce(∪, [simplex_boundary(x) for x in cubes])
    cube_bound_bound = reduce(∪, [simplex_boundary(x) for x in cube_bound])
    cube_bound_bound_bound = reduce(∪, [simplex_boundary(x) for x in cube_bound_bound])
    square_bound = reduce(∪, [simplex_boundary(x) for x in squares])
    square_bound_bound = reduce(∪, [simplex_boundary(x) for x in square_bound])
    edge_bound = reduce(∪, [simplex_boundary(x) for x in edges])
    vertices = setdiff(vertices, cube_bound_bound_bound ∪ square_bound_bound ∪ edge_bound)
    edges = setdiff(edges, cube_bound_bound ∪ square_bound)
    squares = setdiff(squares, cube_bound)

    return (vertices, edges, squares, cubes)
end

# Create a dictionary of all simplices grouped by dimension for a simplicial complex
# represented by its maximal simplices. I dont think we need this
# function dim_dict_scx(max_scx)
    
# end

# This function returns the n-th boundary matrix Dn with coefficients Z/pZ or R of a simplicial complex given as a dimension_dict.
# (p should be a prime. For convenience, p = 0 means R.)
function boundary_matrix(scx, n)
    # scx = [[y for y in x] for x in scx]
    # cases n = 0 and n = dim+1 should be treated separately
    if n == 0
        Dn = zeros(Float64, 1, length(scx[n+1]))
    elseif n == length(scx)
        Dn = zeros(Float64, length(scx[n]), 1)
    else
        Dn = zeros(Float64, length(scx[n]), length(scx[n+1]))
        # otherwise, just replace the zeroes with 1's or -1's at appropriate places
        for j in eachindex(scx[n + 1])
            boundary = collect(simplex_boundary(scx[n + 1][j]))
            for i in eachindex(boundary)
                # println(boundary[i])
                # println(findfirst(==(boundary[i]), scx[n]))
                # if findfirst(==(boundary[i]), scx[n]) === nothing
                #     println(scx[n])
                # end
                Dn[findfirst(==(boundary[i]), scx[n]), j] = 1
            end
        end
    end
    return Dn
end

# determine the mod p Betti numbers of a simplicial complex given via its maximal simplices
function betti_numbers(scx)
    println(scx)
    scx = [[y for y in x] for x in scx]
    scx[2:end] = [[Set(y) for y in x] for x in scx[2:end]]
    # determine the dimension of max_scx
    d = length(scx)
    # initialize the vector of Betti numbers
    betti = []
    # the current (first) boundary matrix
    Dk1 = boundary_matrix(scx, 0)
    # determine the Betti numbers from two consecutive boundary matrices
    for k in 1:d
        # the next boundary matrix
        Dk = boundary_matrix(scx, k)
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