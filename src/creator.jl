function parse_im_data(file_path)
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

function downscale_grid(grid, factor; rule=:majority)
    for f in 1:factor
        dims = size(grid)
        if any(d % 2 != 0 for d in dims)
            error("Grid dimensions must be divisible by 2 at every step.")
        end
        
        # New dimensions after merging 2x2x2 blocks
        new_dims = div.(dims, 2)
        new_grid = zeros(Int, new_dims...)
        
        # Iterate over the new grid positions
        for i in 1:new_dims[1], j in 1:new_dims[2], k in 1:new_dims[3]
            # Extract the 2x2x2 block
            block = grid[
                (i-1)*2+1:i*2,
                (j-1)*2+1:j*2,
                (k-1)*2+1:k*2
            ]
            
            # Apply the rule to determine the new value
            if rule == :majority
                new_grid[i, j, k] = sum(block) > 4 ? 1 : 0
            elseif rule == :max
                new_grid[i, j, k] = maximum(block)
            elseif rule == :min
                new_grid[i, j, k] = minimum(block)
            elseif rule == :random
                new_grid[i, j, k] = block[rand(1:end)]
            else
                error("Unknown rule: $rule")
            end
        end
        
        # Replace the original grid with the new downscaled grid
        grid = new_grid
    end
    
    return grid
end


function cubical_complex(input_grid, factor=0; rule=:min)
    #scale down by a factor
    input_grid = downscale_grid(input_grid, factor, rule=rule)
    vertices = []
    edges = []
    squares = []
    cubes = []

    # Iterate over the 3D grid to extract vertices
    for I in CartesianIndices(input_grid)
        i, j, k = Tuple(I)
        if input_grid[i, j, k] > 0
            # Add the vertex
            push!(vertices, ((i-1)*2, (j-1)*2, (k-1)*2))
        end
    end

    # Extract edges
    for (x,y,z) in vertices
        if (x+2, y, z) in vertices
            push!(edges, (x+1, y, z))
        end
        if (x, Int(y)+2, z) in vertices
            push!(edges, (x, y+1, z))
        end
        if (x, y, z+2) in vertices
            push!(edges, (x, y, z+1))
        end
    end

    # extract squares
    for (x,y,z) in edges
        if isodd(x)
            # check for square in y direction
            if (x-1, y+1, z) in edges && (x+1, y+1, z) in edges && (x, y+2, z) in edges
                push!(squares, (x, y+1, z))
            end
            # check for square in z direction
            if (x-1, y, z+1) in edges && (x+1, y, z+1) in edges && (x, y, z+2) in edges
                push!(squares, (x, y, z+1))
            end
        elseif isodd(y)
            # check for square in x direction
            if (x, y-1, z+1) in edges && (x, y+1, z+1) in edges && (x, y, z+2) in edges
                push!(squares, (x, y, z+1))
            end
        end
    end

    # extract cubes
    for (x,y,z) in squares
        if isodd(x) && isodd(y)
            if (x,y,z+2) in squares && (x,y-1,z+1) in squares && (x,y+1,z+1) in squares && (x-1,y,z+1) in squares && (x+1,y,z+1) in squares
                push!(cubes, (x,y,z+1))
            end
        end
    end

    return (vertices, edges, squares, cubes)
end