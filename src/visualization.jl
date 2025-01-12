using Plots
using LinearAlgebra

function visualize_qcx((vertices, edges, squares, cubes); 
                       show_vertices=true, show_edges=true, 
                       show_squares=false, show_cubes=false)
    # Create an empty 3D plot
    p = plot(legend=false, size=(800, 600), xlabel="X", ylabel="Y", zlabel="Z", camera=(30,30))
    
    # Plot vertices
    if show_vertices
        xs, ys, zs = [v[1]/2 for v in vertices], [v[2]/2 for v in vertices], [v[3]/2 for v in vertices]
        scatter!(p, xs, ys, zs, marker=:circle, color=:red, label="Vertices", markersize=4)
    end
    
    # Plot edges
    if show_edges
        for (x, y, z) in edges
            # Each edge can be drawn as a line between two connected vertices
            if isodd(x)
                plot3d!(p, [(x-1)/2, (x-1)/2+1], [y/2, y/2], [z/2, z/2], linewidth=2, color=:blue, label="Edges")
            end
            if isodd(y)
                plot3d!(p, [x/2, x/2], [(y-1)/2, (y-1)/2+1], [z/2, z/2], linewidth=2, color=:blue, label="Edges")
            end
            if isodd(z)
                plot3d!(p, [x/2, x/2], [y/2, y/2], [(z-1)/2, (z-1)/2+1], linewidth=2, color=:blue, label="Edges")
            end
        end
    end

    # Plot squares
    if show_squares || show_cubes
        for (x, y, z) in squares
            corners_x = []
            corners_y = []
            corners_z = []
            if iseven(x)
                append!(corners_x, [Int(x/2), Int(x/2),Int(x/2),Int(x/2)])
                append!(corners_y, [y, y-1,y-1,y])
                append!(corners_z, [z, z-1,z,z-1])
            elseif iseven(y)
                append!(corners_x, [x, x-1,x-1,x])
                append!(corners_y, [Int(y/2), Int(y/2), Int(y/2), Int(y/2)])
                append!(corners_z, [z, z-1,z,z-1])
            elseif iseven(z)
                append!(corners_x, [x, x-1,x-1,x])
                append!(corners_y, [y, y-1,y,y-1])
                append!(corners_z, [Int(z/2), Int(z/2),Int(z/2),Int(z/2)])
            end
            
            if show_squares
                mesh3d!(p, corners_x, corners_y, corners_z, color=:lightblue, alpha=0.5, label=false, connections= [(1, 3, 2, 4)] )
            end

            if show_cubes
                if iseven(x)
                    if((x+1, y, z) in cubes) || (x-1, y, z) in cubes
                        mesh3d!(p, corners_x, corners_y, corners_z, color=:green, alpha=0.4, label=false, connections= [(1, 3, 2, 4)] )
                    end
                elseif iseven(y)
                    if((x, y+1, z) in cubes) || (x, y-1, z) in cubes
                        mesh3d!(p, corners_x, corners_y, corners_z, color=:green, alpha=0.4, label=false, connections= [(1, 3, 2, 4)] )
                    end
                elseif iseven(z)
                    if((x, y, z+1) in cubes) || (x, y, z-1) in cubes
                        mesh3d!(p, corners_x, corners_y, corners_z, color=:green, alpha=0.4, label=false, connections= [(1, 3, 2, 4)] )
                    end
                end
            end
        end
    end
    
    # Show the plot
    display(p)
end

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
