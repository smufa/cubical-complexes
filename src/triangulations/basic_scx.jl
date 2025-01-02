# Compute the boundary of a simplex sx and represent it as a set of maximal simplices, 
# i.e., just determine the codimenstion 1 faces of sx. 
# E.g. for sx = (1, 2, 3) we get Set[(1, 2), (1, 3), (2, 3)].
function simplex_boundary(sx)
    # determine the codimension 1 simplices in the boundary of the simplex sx
    # 'remove' one entry of a tuple at each specific index by splicing the corresponding subtuples
    cd1faces = Set([(sx[1:k-1]..., sx[k+1:end]...) for k in eachindex(sx)])
    return cd1faces
end

# Create a dictionary of all simplices grouped by dimension for a simplicial complex
# represented by its maximal simplices.
function dim_dict_scx(max_scx)
    # determine the dimension+1 of the complex
    n = maximum([length(sx) for sx in max_scx])
    # set up the dictionary (empty sets of k-tuples in each dimension)
    dimension_dict = Dict((k-1) => Set{NTuple{k, Int}}() for k in 1:n)
    # add all simplices in max_scx to the dictionary
    for sx in max_scx
        push!(dimension_dict[length(sx)-1], sx)
    end
    # fill the dictionary with all simplices from top down
    for k in (n-1):-1:1
        for sx in dimension_dict[k]
            union!(dimension_dict[k-1], simplex_boundary(sx))
        end
    end
    return dimension_dict
end

# find the Euler characteristic of a simplicial complex given as a list of (maximal) simplices (tuples)
function euler_characteristic(max_scx)
    # initialize the Euler characteristic of an empty abstract simplicial complex
    chi = 0
    # build the dictionary
    scx = dim_dict_scx(max_scx)
    # 'sign-count' the number of simplices
    for k in keys(scx)
        chi += iseven(k) ? length(scx[k]) : -length(scx[k])
    end
    return chi
end

# Build a dictionary representing the inverse poset of a simplicial complex given its maximal simplices.
# Each key => value pair is of the form simplex => 'set of codimension one faces'.
function inv_poset_dict_scx(max_scx)
    # initialize the dictionary
    inv_poset_scx = Dict()
    # start with the 'dimension dictionary'
    dim_scx = dim_dict_scx(max_scx)
    # for each dimension >= 0 add the corresponding faces and new keys
    n = maximum(keys(dim_scx)) # dimension of our simplicial complex
    for k in 0:n
        for simplex in dim_scx[k]
            # add the 'k-dimensional' keys
            push!(inv_poset_scx, simplex => Set{NTuple{k, Int}}())
            # dimension 0 gets special treatment
            if k == 0
                continue
            end
            # add simplex => 'set of faces' to the dictionary
            union!(inv_poset_scx[simplex], simplex_boundary(simplex))
        end
    end
    return inv_poset_scx
end

# Build a dictionary representing the poset of a simplicial complex given its maximal simplices.
# Each key => value pair is of the form simplex => 'set of codimension one cofaces'.
function poset_dict_scx(max_scx)
    # initialize the dictionary
    poset_scx = Dict()
    # start with the 'dimension dictionary'
    dim_scx = dim_dict_scx(max_scx)
    # for each dimension >= 0 add the corresponding cofaces and new keys
    n = maximum(keys(dim_scx)) # dimension of our simplicial complex
    for k in 0:n
        for simplex in dim_scx[k]
            # add the 'k-dimensional' keys
            push!(poset_scx, simplex => Set{NTuple{k+2, Int}}())
            # dimension 0 gets special treatment
            if k == 0
                continue
            end
            # add face => 'set of simplices' to the dictionary
            for bsx in simplex_boundary(simplex)
                push!(poset_scx[bsx], simplex)
            end
        end
    end
    return poset_scx
end

# Determine the open star of a simplex in a simplicial complex given as a poset dictionary.
function star(simplex, poset_scx)
    # initialize the star
    star = Set{Tuple}()
    # add the first simplex to it
    push!(star, simplex)
    # and add all of the faces of that simplex
	# initialize the queue as the set of codimenstion 1 faces
    queue = poset_scx[simplex]
    while !isempty(queue)
		# add all faces from the current queue
        union!(star, queue)
        newqueue = Set{Tuple}()
		# make a 'new queue' with codimension 1 faces of simplices in the current queue
        for sx in queue
            union!(newqueue, poset_scx[sx])
        end
		# replace the queue with the 'new queue'
        queue = newqueue
    end
    return star
end

# Determine the link of a simplex in a simplicial complex given as a poset dictionary.
function link(simplex, poset_scx)
    # determine the (open) star of the simplex
    link = star(simplex, poset_scx)
    # remove the simplex from the open star
    delete!(link, simplex)
    # remove the simplex from each element (tuple) of the open star
    link = Set{Tuple}([filter(ver -> !(ver in simplex), sx) for sx in link])
    return link
end