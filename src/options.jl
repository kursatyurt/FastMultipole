struct Options{TF,TI1,TI2}
    expansion_order::TI1
    n_per_branch::TI2 # max number of elements in a branch
    theta::TF # cutoff radius squared for calculating influence directly
    scale_cell_radius::TF # scales the radius of each cell
    scale_body_radius::TF # scales the radius of each body
end

function Options(; expansion_order=4, n_per_branch=50, theta=4.0, scale_cell_radius=1.00001, scale_body_radius=1.0)
    return Options(expansion_order, n_per_branch, theta, scale_cell_radius, scale_body_radius)
end