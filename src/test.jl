include("creator.jl")
include("geometric_objects.jl")
include("simplicial_homology_q.jl")
using Base.Iterators: filter

function get_complexes(initial_grid, steps)
    complexes = []
    for threshold in 1:steps
        push!(complexes, cubical_complex(initial_grid, threshold))
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

using Plots

function plot_persistent_homology(ph_data; infinite_value=1e3)
    # Create arrays to store birth and death points for each dimension
    all_dims = keys(ph_data)
    colors = [:blue, :green, :red, :orange, :purple]  # Color for each dimension
    max_dim = maximum(all_dims)
    
    plt = plot(title="Persistence Diagram", xlabel="Birth", ylabel="Death", legend=:topright)

    for (dim, points) in ph_data
        if isempty(points)
            continue
        end

        births, deaths = [], []
        for (birth, death) in points
            push!(births, birth)
            if death == Inf
                push!(deaths, infinite_value)  # Replace Inf with a large value
            else
                push!(deaths, death)
            end
        end

        # Scatter plot for current dimension
        scatter!(plt, births, deaths, label="Dim $dim", color=colors[dim+1])
    end

    # Add diagonal line to indicate birth = death
    plot!(plt, [0, infinite_value], [0, infinite_value], lw=2, linestyle=:dash, label="Diagonal", color=:black)

    return plt
end


torus = create_torus_grid_3D(5)

(vertices, edges, squares, cubes) = cubical_complex(torus)
println((vertices, edges, squares, cubes))