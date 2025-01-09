include("geometric_objects.jl")
include("simplicial_homology_q.jl")
using Base.Iterators
using Plots

function parse_im_data(file_path, scale_factor=4)
    # Open the file for reading
    open(file_path, "r") do fid
        # Skip the 1024-byte header
        seek(fid, 1024)

        # Read the binary data from the file
        hex_array = read(fid)

        # Convert the binary data into an array of integers
        integer_array = Int.(hex_array)

        # Reshape the array into a 3D grid of dimensions (128, 128, 128)
        data3d = reshape(integer_array, (128, 128, 128))

        return data3d
    end
end

function cubical_complex(input_grid, inclusion_threshold=0)
    vertices = []
    edges = []
    squares = []
    cubes = []

    # Iterate over the 3D grid to extract vertices
    for I in CartesianIndices(input_grid)
        i, j, k = Tuple(I)
        if input_grid[i, j, k] > inclusion_threshold
            # Add the vertex
            push!(vertices, (i, j, k))
        end
    end

    # Extract edges, squares, and cubes
    for I in CartesianIndices(input_grid)
        i, j, k = Tuple(I)
        if input_grid[i, j, k] > inclusion_threshold
            # Check for edges along each axis
            if i < size(input_grid, 1) && input_grid[i + 1, j, k] > inclusion_threshold
                push!(edges, ((i, j, k), (i + 1, j, k)))
            end
            if j < size(input_grid, 2) && input_grid[i, j + 1, k] > inclusion_threshold
                push!(edges, ((i, j, k), (i, j + 1, k)))
            end
            if k < size(input_grid, 3) && input_grid[i, j, k + 1] > inclusion_threshold
                push!(edges, ((i, j, k), (i, j, k + 1)))
            end

            # Check for squares in the XY, YZ, and XZ planes
            if i < size(input_grid, 1) && j < size(input_grid, 2) &&
               input_grid[i + 1, j, k] > inclusion_threshold && input_grid[i, j + 1, k] > inclusion_threshold && input_grid[i + 1, j + 1, k] > inclusion_threshold
                push!(squares, ((i, j, k), (i + 1, j, k), (i, j + 1, k), (i + 1, j + 1, k)))
            end
            if j < size(input_grid, 2) && k < size(input_grid, 3) &&
               input_grid[i, j + 1, k] > inclusion_threshold && input_grid[i, j, k + 1] > inclusion_threshold && input_grid[i, j + 1, k + 1] > inclusion_threshold
                push!(squares, ((i, j, k), (i, j + 1, k), (i, j, k + 1), (i, j + 1, k + 1)))
            end
            if i < size(input_grid, 1) && k < size(input_grid, 3) &&
               input_grid[i + 1, j, k] > inclusion_threshold && input_grid[i, j, k + 1] > inclusion_threshold && input_grid[i + 1, j, k + 1] > inclusion_threshold
                push!(squares, ((i, j, k), (i + 1, j, k), (i, j, k + 1), (i + 1, j, k + 1)))
            end

            # Check for cubes
            if i < size(input_grid, 1) && j < size(input_grid, 2) && k < size(input_grid, 3) &&
               input_grid[i + 1, j, k] > inclusion_threshold && input_grid[i, j + 1, k] > inclusion_threshold && input_grid[i, j, k + 1] > inclusion_threshold &&
               input_grid[i + 1, j + 1, k] > inclusion_threshold && input_grid[i + 1, j, k + 1] > inclusion_threshold && input_grid[i, j + 1, k + 1] > inclusion_threshold &&
               input_grid[i + 1, j + 1, k + 1] > inclusion_threshold
                push!(cubes, ((i, j, k), (i + 1, j, k), (i, j + 1, k), (i, j, k + 1),
                              (i + 1, j + 1, k), (i + 1, j, k + 1), (i, j + 1, k + 1), (i + 1, j + 1, k + 1)))
            end
        end
    end

    return (vertices, edges, squares, cubes)
end

function visualize_cubical_complex(vertices, edges; show_edges=true, edge_color=:blue)
    # Extract x, y, z coordinates for vertices
    x = [v[1] for v in vertices]
    y = [v[2] for v in vertices]
    z = [v[3] for v in vertices]

    # Create a scatter plot for vertices
    scatter3d(x, y, z, markersize=5, label="Vertices", color=:red, legend=:top)

    if show_edges
        # Plot edges
        for e in edges
            x_edge = [e[1][1], e[2][1]]
            y_edge = [e[1][2], e[2][2]]
            z_edge = [e[1][3], e[2][3]]
            plot3d!(x_edge, y_edge, z_edge, color=edge_color, linewidth=2, label=false)
        end
    end

    display(plot!())
end
