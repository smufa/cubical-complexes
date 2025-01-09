include("cubical_complex.jl")
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

initial_grid = create_torus_grid_3D(10)
complexes = get_complexes(initial_grid, 10)

for complex in complexes
    (vertices, edges, squares, cubes) = complex

    visualize_cubical_complex(vertices, edges)

    println("vertices: ", length(vertices))
    println("edges: ", length(edges))
    println("squares: ", length(squares))
    println("cubes: ", length(cubes))

    println("")
end