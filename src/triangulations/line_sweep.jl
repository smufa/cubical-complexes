
using Plots

function do_segments_intersect(AB, CD)
    # form the matrix of the corresponding linear system
    A = [AB[:, 2]-AB[:, 1] CD[:, 1]-CD[:, 2]]
    # form the right-hand side
    b = CD[:, 1]-AB[:, 1]
    # solve the linear system
    t = A\b
    return t[1]>=0 && t[1]<=1 && t[2]>=0 && t[2]<=1
end

function line_sweep(points)
    # Order the points wrt. to the x-coordinate and return the corresponding permutation.
    # (Our function will return an abstract simplicial complex.)
    sorted = sortperm(points, by = x -> x[1])
    # first three points determine the first triangle (first triangle is needed to start the algorithm)
    abstract_triangulation = [(sorted[1],), (sorted[2],), (sorted[1], sorted[2]), (sorted[3],), (sorted[1], sorted[3]), (sorted[2], sorted[3]), (sorted[1], sorted[2], sorted[3])]
    boundary_vertices = [sorted[1], sorted[2], sorted[3]]
    # and we start the sweep
    for vertex in sorted[4:end]
        # find boundary edges that intersect the segment from vertex to any of the boudary vertices
        n = length(boundary_vertices)
        # store indices of visible vertices and visible edges in a boolean list
        visible_vert_ind = trues(n)
        visible_edge_ind = trues(n)
        for bv_ind in 1:n # loop through boundary vertices (bvertex)
            for k in 1:n # loop through boundary edges (bedge)
                # segment from bvertex to next boundary vertex
                CD = [collect(points[boundary_vertices[k]]) collect(points[boundary_vertices[k != n ? k+1 : 1]])]
                # skip bedges, for which bvertex is one of the vertices, from vertex-bvertex intersections
                if !(bv_ind == k)
                    # segment from vertex to the midpoint of the bedge
                    AB = [collect(points[vertex]) (collect(points[boundary_vertices[bv_ind]]) + collect(points[boundary_vertices[bv_ind != n ? bv_ind+1 : 1]]))/2]
                    # remove bedge bv_ind from visible edges if AB and CD intersect
                    # (&= is used since it should remain removed if it was removed in an earlier iteration)
                    visible_edge_ind[bv_ind] &= !do_segments_intersect(AB, CD)
                end
                # skip bvertices, for which bvertex is one of the vertices, from vertex-bvertex intersections
                if !(bv_ind == k || bv_ind == (k != n ? k+1 : 1))
                    # segment from vertex to bvertex
                    AB = [collect(points[vertex]) collect(points[boundary_vertices[bv_ind]])]
                    # remove bvertex bv_ind from visible vertices if AB and CD intersect
                    visible_vert_ind[bv_ind] &= !do_segments_intersect(AB, CD)
                end
            end
        end
        # add the current vertex to abstract triangulation
        push!(abstract_triangulation, (vertex,))
        # add all edges from vertex to visible vertices and all triangles from vertex 'to' visible edges to the abstract triangulation
        for k in 1:n
            if visible_vert_ind[k]
                push!(abstract_triangulation, (boundary_vertices[k]..., vertex))
            end
            if visible_edge_ind[k]
                push!(abstract_triangulation, (boundary_vertices[[k (k != n ? k+1 : 1)]]..., vertex))
            end
        end

        # The visible vertices (apart from the 'first' and the 'last' one) have just become invisible.
        # Since the lists of boundary vertices and boundary edges are circular lists
        # we deem that the 'first' and the 'last' vertex are the trues adjacent to a false 
        # in boundary vertices list. If there are no false-s in the boundary vertices list, 
        # there is just one false in the boundary edges list and the 'first' and the 'last' 
        # vertex are the ones adjacent to that false.

        # Let's find the true-false and the false-true occurences.
        # (This will be the index of a true adjacent to a false.)
        tf_ind = nothing # true-to-false index, the 'last' one
        ft_ind = nothing # false-to-true index, the 'first' one
        current = visible_vert_ind[1]
        for k in 2:n
            if current != visible_vert_ind[k]
                current ? tf_ind = k-1 : ft_ind = k
            end
            current = visible_vert_ind[k]
        end
        if isnothing(tf_ind)
            if isnothing(ft_ind)
                only_false = findfirst(x -> !x, visible_edge_ind)
                tf_ind = only_false
                ft_ind = (only_false != n ? only_false+1 : 1)
            else
                tf_ind = n
            end
        elseif isnothing(ft_ind)
            ft_ind = 1
        end
        # delete the vertices, which have just become invisible from the boundary vertices list
        visible_vert_ind[[ft_ind tf_ind]] = [false false]
        deleteat!(boundary_vertices, visible_vert_ind)
        # insert the current vertex to boundary vertices list at the correct position
        if ft_ind < tf_ind
            insert!(boundary_vertices, ft_ind+1, vertex)
        elseif tf_ind == 1
            push!(boundary_vertices, vertex)
        else
            insert!(boundary_vertices, tf_ind-1, vertex)
        end
    end
    return abstract_triangulation
end

function visual_test_ls(num_points)
    # generate uniformly random points in the square [0,1]x[0,1]
    points = [(rand(), rand()) for n in 1:num_points]
    # perform a line sweep trianulation
    triangulation = line_sweep(points)
    # clear the plot
    plot()
    # draw the vertices and the edges
    for simplex in triangulation
        if length(simplex) == 1
            # plot a vertex
            scatter!(points[simplex[1]], markersize=2.5, label="")
        elseif length(simplex) == 2
            # plot an edge
            plot!([points[simplex[1]][1], points[simplex[2]][1]], [points[simplex[1]][2], points[simplex[2]][2]], label="")
        end
    end
    # show the plot
    display(plot!())
end

visual_test_ls(20)
