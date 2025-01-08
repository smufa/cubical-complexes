function create_torus_grid_3D(n)
    torus = zeros(Int, n, n, n)  # Create a 3D grid of size n x n x n

    for i in 1:n
        for j in 1:n
            for k in 1:n
                if (i == 1 || j == 1 || i == n || j == n
                    || i == 2 || j == 2 || i == n-1 || j == n-1)
                    torus[i, j, 1] = 1
                    torus[i, j, 2] = 1
                end
            end
        end
    end

    return torus
end

function create_torus_grid_2D(n)
    torus = zeros(Int, n, n, 1)  # Create a 3D grid of size n x n x 1

    for i in 1:n
        for j in 1:n
            for k in 1:n
                if (i == 1 || j == 1 || i == n || j == n
                    || i == 2 || j == 2 || i == n-1 || j == n-1)
                    torus[i, j, 1] = 1
                end
            end
        end
    end

    return torus
end

function create_klein_bottle_boundary(n)
    edges = []
    squares = []

    # Add edges
    for i in 1:n, j in 1:n
        # Horizontal edges (wrap horizontally)
        horizontal_edge = ((i, j), (i % n + 1, j))
        push!(edges, horizontal_edge)

        # Vertical edges (wrap vertically)
        vertical_edge = ((i, j), (i, j % n + 1))
        push!(edges, vertical_edge)
    end

    # Add squares
    for i in 1:n, j in 1:n
        # Define corners of the square
        v1 = (i, j)
        v2 = (i % n + 1, j)
        v3 = (i, j % n + 1)
        v4 = (i % n + 1, j % n + 1)

        # Add the square, represented by its 4 edges
        square = (v1, v2, v3, v4)
        push!(squares, square)
    end

    return (edges, squares)
end

function create_projective_plane_boundary(n)
    edges = []
    squares = []

    # Add edges
    for i in 1:n, j in 1:n
        # Horizontal edges with a twist
        horizontal_edge = ((i, j), (n - i + 1, j))  # Opposite edges reversed
        push!(edges, horizontal_edge)

        # Vertical edges with a twist
        vertical_edge = ((i, j), (i, n - j + 1))  # Opposite edges reversed
        push!(edges, vertical_edge)
    end

    # Add squares
    for i in 1:n, j in 1:n
        # Define corners of the square
        v1 = (i, j)
        v2 = (n - i + 1, j)  # Opposite vertex with horizontal twist
        v3 = (i, n - j + 1)  # Opposite vertex with vertical twist
        v4 = (n - i + 1, n - j + 1)  # Diagonal twist

        # Add square
        push!(squares, (v1, v2, v3, v4))
    end

    return (edges, squares)
end
