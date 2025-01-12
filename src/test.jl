include("creator.jl")
include("geometric_objects.jl")
include("simplicial_homology_q.jl")
include("visualization.jl")
using Base.Iterators: filter

function get_complexes(initial_grid, steps)
    complexes = []
    for threshold in 0:steps
        push!(complexes, cubical_complex(initial_grid, threshold, rule=:min))
    end

    return complexes
end

function persistent_homology(filtration)
    # filtration: A list of simplicial complexes at different levels of the filtration
    # Each simplicial complex is represented as a tuple of (vertices, edges, squares, cubes, etc.)
    
    persistence_pairs = Dict()  # Store birth and death pairs for each dimension
    for dim in 0:(length(filtration[1]) - 1)  # For each dimension of homology
        persistence_pairs[dim] = []
    end
    
    # Track active homology classes using a dictionary
    active_classes = Dict{Int, Dict}()
    
    for dim in 0:(length(filtration[1]) - 1)
        active_classes[dim] = Dict()
    end
    
    # Iterate over filtration levels
    for i in 0:length(filtration)-1
        scx = filtration[i+1]
        # Compute the boundary matrix for each dimension
        for dim in 0:(length(scx) - 1)
            D = boundary_matrix(scx, dim)
            
            # Perform Gaussian elimination to track births and deaths
            rows, cols = size(D)
            for col in 1:cols
                pivot_row = findfirst(row -> D[row, col] != 0, 1:rows)
                if pivot_row !== nothing
                    # This column has a pivot => death of a homology class
                    if haskey(active_classes[dim], pivot_row)
                        birth = active_classes[dim][pivot_row]
                        push!(persistence_pairs[dim], (birth, i))
                        delete!(active_classes[dim], pivot_row)
                    end
                else
                    # No pivot => birth of a new homology class
                    active_classes[dim][col] = i
                end
            end
        end
    end
    
    # Any remaining active classes represent homology that persists to infinity
    for dim in 0:(length(filtration[1]) - 1)
        for (_, birth) in active_classes[dim]
            push!(persistence_pairs[dim], (birth, Inf))
        end
    end
    
    return persistence_pairs
end


# Now the more intersting object:
intersting_object = create_hollow_cube_with_holes_and_hollow_inner_cube(16, thickness=1, inner_edge_length=4, hole_edge_length=1)

cx = cubical_complex(intersting_object, 1; rule=:max)

println(betti_numbers(cx))

visualize_qcx(cx, show_squares=false)