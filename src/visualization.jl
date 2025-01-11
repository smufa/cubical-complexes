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

function visualize_qcx((vertices, edges, squares, cubes); show_vertices=true, show_edges=true, show_cubes=false, show_squares=false, edge_color=:blue, square_color=:lightblue, cube_color=:lightgreen)
    p = plot!()
    verticies_x = []
    verticies_y = []
    verticies_z = []
    if show_vertices
        for (x, y, z) in vertices
            push!(verticies_x, (x/2))
            push!(verticies_y, (y/2))
            push!(verticies_z, (z/2))
        end
        scatter3d(p, verticies_x, verticies_y, verticies_z, markersize=5, label="Vertices", color=:red)
    end

    if show_edges
        for (x, y, z) in edges
            if isodd(x)
                x_edge = (Int(((x-1)/2)), Int(((x-1)/2) + 1))
                y_edge = (Int((y/2)), Int((y/2)))
                z_edge = (Int((z/2)), Int((z/2)))
            elseif isodd(y)
                x_edge = (Int((x/2)), Int((x/2)))
                y_edge = (Int(((y - 1)/2)), Int(((y - 1)/2) + 1))
                z_edge = (Int((z/2)), Int((z/2)))
            else
                x_edge = (Int((x/2)), Int((x/2)))
                y_edge = (Int((y/2)), Int((y/2)))
                z_edge = (Int(((z - 1)/2)), Int(((z - 1)/2) + 1))
            end
            plot3d!(p, x_edge, y_edge, z_edge, color=edge_color, linewidth=2, label=false)
        end
    end

    if show_squares

    end

    if show_cubes

    end

    display(p)
end
