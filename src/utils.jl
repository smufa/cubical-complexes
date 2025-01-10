using Interpolations

function resample_grid(grid::Array{T, 3}, factor::Int) where T
    size_x, size_y, size_z = size(grid)
    new_size = div.(size(grid), factor)
    
    # Create interpolation object
    itp = interpolate(grid, BSpline(Linear()), OnGrid())
    scale_x = range(1, stop=size_x, length=new_size[1])
    scale_y = range(1, stop=size_y, length=new_size[2])
    scale_z = range(1, stop=size_z, length=new_size[3])
    
    # Resample
    return [itp[x, y, z] for x in scale_x, y in scale_y, z in scale_z]
end