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

function cubical_complex(input_grid)
    verticies_to_check = []
    edges_to_check = []
    squares_to_check = []

    vertices = []
    edges = []
    squares = []
    cubes = []

    grid_size = 0

    # Iterate over the 3D grid to extract vertices
    for I in CartesianIndices(input_grid)
        i, j, k = Tuple(I)
        if input_grid[i, j, k] > 0
            # Add the vertex
            push!(vertices, (i*2, j*2, k*2))
            push!(verticies_to_check, (i, j, k))
        end
        grid_size = i
    end

    # Extract edges
    for (i, j, k) in verticies_to_check
        if i == grid_size || j == grid_size || k == grid_size
            break
        end
        if input_grid[i + 1, j, k] > 0
            push!(edges, (i*2 + 1, j*2, k*2))
            push!(edges_to_check, (i, j, k))
        elseif input_grid[i, j + 1, k] > 0
            push!(edges, (i*2, j*2 + 1, k*2))
            push!(edges_to_check, (i, j, k))
        else  input_grid[i, j, k + 1] > 0
            push!(edges, (i*2, j*2, k*2 + 1))
            push!(edges_to_check, (i, j, k))
        end
    end

    # extract squares
    for (i, j, k) in edges_to_check
        if i == grid_size || j == grid_size || k == grid_size
            break
        end
        if input_grid[i + 1, j + 1, k] > 0 && input_grid[i + 1, j, k] > 0 && input_grid[i, j + 1, k] > 0
            push!(squares, (i*2 + 1, j*2 + 1, k*2))
            push!(squares_to_check, (i, j, k))
        elseif  input_grid[i + 1, j, k + 1] > 0 && input_grid[i + 1, j, k] > 0 && input_grid[i, j, k + 1] > 0
            push!(squares, (i*2 + 1, j*2, k*2 + 1))
            push!(squares_to_check, (i, j, k))
        else  input_grid[i, j + 1, k + 1] > 0 && input_grid[i, j + 1, k] > 0 && input_grid[i, j, k + 1] > 0
            push!(squares, (i*2, j*2 + 1, k*2 + 1))
            push!(squares_to_check, (i, j, k))
        end
    end

    # extract cubes
    for (i, j, k) in edges_to_check
        if input_grid[i + 1, j + 1, k + 1] > 0 && input_grid[i + 1, j, k] > 0 && input_grid[i, j + 1, k] > 0 && input_grid[i, j, k + 1] > 0 && input_grid[i, j + 1, k + 1] > 0 && input_grid[i+ 1, j, k + 1] > 0 && input_grid[i + 1, j + 1, k] > 0
            push!(cubes, (i*2 + 1, j*2 + 1, k*2 + 1))
        end
    end

    return (vertices, edges, squares, cubes)
end