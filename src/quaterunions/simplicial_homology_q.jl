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
function boundary_matrix(dimension_dict, n, p = 0)
    # # check if p is 0 or prime
    # if p != 0 && !isprime(p)
    #     throw(DomainError(p, "p should be 0 or a prime."))
    # end
    # # cases n = 0 and n = dim+1 should be treated separately
    # if n == 0
    #     Dn = zeros(p==0 ? Float64 : Mod{p}, 1, length(dimension_dict[n]))
    # elseif n == maximum(keys(dimension_dict))+1
    #     Dn = zeros(p==0 ? Float64 : Mod{p}, length(dimension_dict[n-1]), 1)
    # else
    #     Dn = zeros(p==0 ? Float64 : Mod{p}, length(dimension_dict[n-1]), length(dimension_dict[n]))
    #     # otherwise, just replace the zeroes with 1's or -1's at appropriate places
    #     for j in eachindex(dimension_dict[n])
    #         boundary = simplex_boundary(dimension_dict[n][j])
    #         for i in eachindex(boundary)
    #             Dn[findfirst(==(boundary[i]), dimension_dict[n-1]), j] = (-1)^i
    #         end
    #     end
    # end
    # return Dn
end

# determine the mod p Betti numbers of a simplicial complex given via its maximal simplices
function betti_numbers(max_scx, p = 0)
    # TODO
end

# This function performs a simultaneous reduction on matrices A=D_n and B=D_{n+1}
# to determine the number of generators of H_n. It implicitly assumes that AB=0.
function simultaneous_reduce(A, B)
    # TODO
end