using LinearAlgebra

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

function create_sphere_grid_3D(n::Int; radius=4)
    grid = zeros(Int, n, n, n)  # Initialize all values to 0
    center = (div(n, 2), div(n, 2), div(n, 2))
    radius = radius > 0 ? radius : div(n, 3)
    for i in 1:n, j in 1:n, k in 1:n
        if norm((i, j, k) .- center) <= radius
            grid[i, j, k] = 1  # Mark as present
        end
    end
    return grid
end

function create_hollow_cube_grid_3D(n::Int; thickness=1)
    grid = ones(Int, n, n, n)  # Initialize all values to 0
    grid[1+thickness:(end-thickness), 1+thickness:(end-thickness), 1+thickness:(end-thickness)] .= 0
    return grid
end

function create_twist_grid_3D(n::Int; twist_factor=5)
    grid = zeros(Int, n, n, n)  # Initialize all values to 0
    for i in 1:n, j in 1:n, k in 1:n
        if abs(i - div(n, 2)) <= twist_factor * sin(2 * π * (j / n)) &&
           abs(k - div(n, 2)) <= twist_factor * sin(2 * π * (j / n))
            grid[i, j, k] = 1  # Mark as present
        end
    end
    return grid
end

function create_random_grid_3D(n::Int; density=0.1)
    grid = zeros(Int, n, n, n)  # Initialize all values to 0
    for i in 1:n, j in 1:n, k in 1:n
        if rand() < density
            grid[i, j, k] = 1  # Mark as present
        end
    end
    return grid
end

function create_hollow_cube_with_holes_and_hollow_inner_cube(
    n::Int; thickness=1, inner_edge_length=3, inner_thickness=1, hole_edge_length=1
)
    # Initialize the grid with all ones (outer cube)
    grid = ones(Int, n, n, n)
    
    # Create the hollow part of the outer cube
    grid[1+thickness:(end-thickness), 1+thickness:(end-thickness), 1+thickness:(end-thickness)] .= 0
    
    # Add holes to each face of the outer cube
    hole_start = div(n - hole_edge_length, 2) + 1
    hole_end = hole_start + hole_edge_length - 1
    
    # Front and back faces (x-y plane at z=1 and z=n)
    grid[:, hole_start:hole_end, 1] .= 0
    grid[:, hole_start:hole_end, n] .= 0
    
    # Left and right faces (y-z plane at x=1 and x=n)
    grid[1, hole_start:hole_end, :] .= 0
    grid[n, hole_start:hole_end, :] .= 0
    
    # Top and bottom faces (x-z plane at y=1 and y=n)
    grid[hole_start:hole_end, 1, :] .= 0
    grid[hole_start:hole_end, n, :] .= 0
    
    # Add the inner hollow cube
    center_start = div(n - inner_edge_length, 2) + 1
    center_end = center_start + inner_edge_length - 1
    grid[center_start:center_end, center_start:center_end, center_start:center_end] .= 1  # Inner cube outer shell
    
    # Hollow out the inner cube
    inner_hollow_start = center_start + inner_thickness
    inner_hollow_end = center_end - inner_thickness
    if inner_hollow_start <= inner_hollow_end
        grid[inner_hollow_start:inner_hollow_end, inner_hollow_start:inner_hollow_end, inner_hollow_start:inner_hollow_end] .= 0
    end
    
    return grid
end
