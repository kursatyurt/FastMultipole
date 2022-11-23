const BRANCH_TYPE = Float64

struct Branch{TF,TFV1<:AbstractVector{TF},TFV3<:AbstractVector{TF},TF2<:NTuple{4,<:AbstractVector{Complex{TF}}},TI8<:Integer,TI32<:Integer,TIV32<:Vector{TI32}}
    n_branches::TI8          # number of child branches
    n_bodies::TIV32          # number of descendent bodies
    first_branch::TI32        # index of the first branch
    first_body::TIV32       # index of the first element
    center::TFV3 # center of the branch (i=1:3) and its radius (i=4)
    radius::TFV1              # side lengths of the rectangle encapsulating the branch
    multipole_expansion::TF2 # multipole expansion coefficients
    local_expansion::TF2     # local expansion coefficients
end

struct Tree{TBV,TI1,TI2}
    branches::TBV        # a vector of `Branch` objects composing the tree
    expansion_order::TI1
    n_per_branch::TI2    # max number of bodies in a leaf
end

"""
    Tree(elements; expansion_order=2, n_per_branch=1)

Constructs an octree of the provided element objects.

# Inputs

- `elements`- a struct containing the following members:

    * `bodies::Array{Float64,2}`- 3+4+mxN array containing element positions, strengths, and m other values that must be sorted into the octree
    * `potential::Array{Float64,2}`- 4+12+36xN array of the potential, Jacobian, and Hessian that are reset every iteration (and don't require sorting)
    * `velocity::Array{Float64,2}`- 3xN array of the velocity vectors at each element, reset every iteration, and calculated in post-processing
    * `direct!::Function`- function calculates the direct influence of the body at the specified location
    * `B2M!::Function`- function converts the body's influence into a multipole expansion
"""
function Tree(elements_tuple::Tuple, options::Options)
    # read options
    expansion_order = options.expansion_order
    n_per_branch = options.n_per_branch
    scale_cell_radius = options.scale_cell_radius
    scale_body_radius = options.scale_body_radius

    # initialize objects
    bodies_list = [elements.bodies for elements in elements_tuple]
    buffer_list = [similar(bodies) for bodies in bodies_list]
    index_list = [elements.index for elements in elements_tuple]
    # branch_type = typeof() # TODO: make typing implicit
    branches = Vector{Branch{BRANCH_TYPE,MVector{1, BRANCH_TYPE},MVector{3, BRANCH_TYPE},NTuple{4,MVector{((expansion_order+1) * (expansion_order+2)) >> 1, Complex{BRANCH_TYPE}}},Int8,Int32,Vector{Int32}}}(undef,1)

    # recursively build branches
    i_start = ones(Int32,length(elements_tuple))
    i_end = [Int32(size(bodies)[2]) for bodies in bodies_list]
    i_branch = 1
    center, radius = center_radius(elements_tuple,(i_start[i]:i_end[i] for i in 1:length(elements_tuple)), scale_body_radius, scale_cell_radius)
    level = 0
    branch!(branches, bodies_list, buffer_list, index_list, i_start, i_end, i_branch, center, radius, level, expansion_order, n_per_branch)

    # assemble tree
    tree = Tree(branches, Int16(expansion_order), Int32(n_per_branch))

    # update branch centers/radii
    # adjust_center_radius!(tree, elements_tuple, options)

    return tree
end

function branch!(branches, bodies_list, buffer_list, index_list, i_start, i_end, i_branch, center, radius, level, expansion_order, n_per_branch)
    n_branches = Int8(0)
    n_bodies = i_end - i_start .+ Int32(1)
    n_types = length(i_start)
    multipole_expansion = (
        MVector{((expansion_order+1) * (expansion_order+2)) >> 1,Complex{typeof(branches).parameters[1].parameters[1]}}(undef),
        MVector{((expansion_order+1) * (expansion_order+2)) >> 1,Complex{typeof(branches).parameters[1].parameters[1]}}(undef),
        MVector{((expansion_order+1) * (expansion_order+2)) >> 1,Complex{typeof(branches).parameters[1].parameters[1]}}(undef),
        MVector{((expansion_order+1) * (expansion_order+2)) >> 1,Complex{typeof(branches).parameters[1].parameters[1]}}(undef)
    )
    local_expansion = (
        MVector{((expansion_order+1) * (expansion_order+2)) >> 1,Complex{typeof(branches).parameters[1].parameters[1]}}(undef),
        MVector{((expansion_order+1) * (expansion_order+2)) >> 1,Complex{typeof(branches).parameters[1].parameters[1]}}(undef),
        MVector{((expansion_order+1) * (expansion_order+2)) >> 1,Complex{typeof(branches).parameters[1].parameters[1]}}(undef),
        MVector{((expansion_order+1) * (expansion_order+2)) >> 1,Complex{typeof(branches).parameters[1].parameters[1]}}(undef)
    )
    if sum(n_bodies) <= n_per_branch # checks for all element structs => branch is a leaf; no new branches needed
        i_child = Int32(-1)
        branch = Branch(n_branches, n_bodies, i_child, i_start, center, radius, multipole_expansion, local_expansion)
        branches[i_branch] = branch
        return nothing
    else # not a leaf; branch children
        # count elements in each octant
        octant_attendance = zeros(Int32, 8, n_types)
        for i_type in 1:n_types # loop over each type
            for i_body in i_start[i_type]:i_end[i_type] # loop over each body
                x = bodies_list[i_type][1:3,i_body]
                i_octant = get_octant(x, center)
                octant_attendance[i_octant, i_type] += Int32(1)
            end
        end

        # determine number of children
        for i_octant in 1:8
            n_branches += true in (octant_attendance[i_octant,:] .> 0) # check for elements of any type
        end
        
        # create branch
        i_child = Int32(length(branches) + 1)
        branches[i_branch] = Branch(n_branches, n_bodies, i_child, i_start, center, radius, multipole_expansion, local_expansion)

        # sort bodies
        ## write offsets
        offsets = zeros(Int32,8,n_types)
        offsets[1,:] .= i_start
        for i=2:8
            offsets[i,:] .= offsets[i-1,:] + octant_attendance[i-1,:]
        end

        ## create counter
        counter = deepcopy(offsets)

        ## sort element indices into the buffer
        for i_type in 1:n_types
            for i_body in i_start[i_type]:i_end[i_type]
                x = bodies_list[i_type][1:3,i_body]
                i_octant = get_octant(x, center)
                buffer_list[i_type][:,counter[i_octant,i_type]] .= bodies_list[i_type][:,i_body]
                index_list[i_type][i_body] = counter[i_octant,i_type]
                counter[i_octant,i_type] += 1
            end
            # place sorted bodies
            bodies_list[i_type][:,i_start[i_type]:i_end[i_type]] .= buffer_list[i_type][:,i_start[i_type]:i_end[i_type]]
        end

        # recursively build new branches
        prev_branches = length(branches)
        resize!(branches, prev_branches + n_branches)

        i_tape = 1
        for i_octant in Int8(0):Int8(7) # check each octant for members
            if true in (octant_attendance[i_octant+Int8(1),:] .> 0) # if octant is populated
                child_i_start = offsets[i_octant+Int8(1),:] # index of first member for all element types
                child_i_end = child_i_start + octant_attendance[i_octant+Int8(1),:] .- Int32(1) # if attendence is 0, the end index will be less than the start
                child_radius = radius / 2#(1 << (level + 1))
                child_center = deepcopy(center)
                for d in Int8(0):Int8(2)
                    child_center[d + Int8(1)] += child_radius[1] * (((i_octant & Int8(1) << d) >> d) * Int8(2) - Int8(1))
                end
                branch!(branches, bodies_list, buffer_list, index_list, child_i_start, child_i_end, i_tape + prev_branches, child_center, child_radius, level + 1, expansion_order, n_per_branch)
                i_tape += 1
            end
        end
    end
end

function center_radius(elements_tuple::Tuple, indices, scale_body_radius, scale_cell_radius)
    x_min = MVector{3}(elements_tuple[1].bodies[i_POSITION,1])
    x_max = MVector{3}(elements_tuple[1].bodies[i_POSITION,1])
    for (i_type,index) in enumerate(indices)
        for i in index
            x = elements_tuple[i_type].bodies[i_POSITION,i]
            r = elements_tuple[i_type].bodies[i_SIGMA,i] * scale_body_radius
            for dim in 1:3
                if x[dim]-r < x_min[dim]
                    x_min[dim] = x[dim]-r
                elseif x[dim]+r > x_max[dim]
                    x_max[dim] = x[dim]+r
                end
            end
        end
    end
    center = (x_max + x_min) ./ 2
    @show scale_cell_radius max(x_max-center...)
    # radius = MVector{1}(max(x_max-center...) * scale_cell_radius) # use the infinity norm
    radius = MVector{1}(norm(x_max-center) * scale_cell_radius) # use the L2 norm
    return center, radius
end

function center_radius(branches, i_branch::Int, scale_cell_radius)
    this_branch = branches[i_branch]
    x_min = MVector{3}(this_branch.center)
    x_max = MVector{3}(this_branch.center)
    for child_branch in branches[this_branch.first_branch:this_branch.first_branch+this_branch.n_branches-1]
        x = child_branch.center
        r = child_branch.radius[1]
        for dim in 1:3
            if x[dim]-r < x_min[dim]
                x_min[dim] = x[dim]-r
            elseif x[dim]+r > x_max[dim]
                x_max[dim] = x[dim]+r
            end
        end
    end
    center = (x_max + x_min) ./ 2
    println("branches")
    @show scale_cell_radius max(x_max-center...)
    # radius = MVector{1}(max(x_max-center...) * scale_cell_radius) # use the infinity norm
    radius = MVector{1}(norm(x_max-center)) # use the L2 norm; don't need to scale the radius again (see children)
    return center, radius
end

function adjust_center_radius!(branches, i_branch, elements_tuple, scale_body_radius, scale_cell_radius)
    this_branch = branches[i_branch]

    # recurse over branches
    for i_child in this_branch.first_branch:this_branch.first_branch+this_branch.n_branches-1
        adjust_center_radius!(branches, i_child, elements_tuple, scale_body_radius, scale_cell_radius)
    end

    # adjust this branch if it is a leaf
    if sum(this_branch.n_bodies) < 0
        indices = (this_branch.first_body[i_type]:this_branch.first_body[i_type]+
            this_branch.n_bodies-1 for i_type in 1:length(elements_tuple))
        center, radius = center_radius(elements_tuple, indices, scale_body_radius, scale_cell_radius)
    else # if it is not a leaf
        center, radius = center_radius(branches, i_branch::Int, scale_cell_radius)
    end
    this_branch.center .= center
    this_branch.radius .= radius
    return nothing
end

function adjust_center_radius!(tree, elements_tuple, options)
    branches = tree.branches
    scale_body_radius = options.scale_body_radius
    scale_cell_radius = options.scale_cell_radius
    adjust_center_radius!(branches, 1, elements_tuple, scale_body_radius, scale_cell_radius)
end

@inline function get_octant(x, center)
    i_octant = (UInt8((x[1] > center[1])) + UInt8(x[2] > center[2]) << 0b1 + UInt8(x[3] > center[3]) << 0b10) + 0b1
end

# function n_terms(expansion_order, dimensions)
#     n = 0
#     for order in 0:expansion_order
#         # leverage stars and bars theorem
#         n += binomial(order + dimensions - 1, dimensions - 1)
#     end
#     return n
# end

# function initialize_expansion(expansion_order)
#     return [zeros(Complex{Float64}, ((expansion_order+1) * (expansion_order+2)) >> 1) for _ in 1:4]
# end

# function change_expansion_order!(tree::Tree, new_order)
#     old_n = length(tree.branches[1].multipole_expansion[1])
#     new_n = n_terms(new_order, 3)
#     old_order = tree.expansion_order
#     for branch in tree.branches
#         if  old_order < new_order
#             for i in 1:new_n - old_n
#                 for dim in 1:4
#                     push!(branch.multipole_expansion[dim],0.0)
#                     push!(branch.local_expansion[dim],0.0)
#                 end
#             end
#         elseif old_order > new_order
#             for i in 1:new_n - old_n
#                 for dim in 1:4
#                     pop!(branch.multipole_expansion[dim])
#                     pop!(branch.local_expansion[dim])
#                 end
#             end
#         end
#     end
#     tree.expansion_order = new_order
# end

function reset_expansions!(tree)
    for branch in tree.branches
        branch.multipole_expansion .*= 0
        branch.local_expansion .*= 0
    end
end

Base.eltype(branch::Branch) = typeof(branch).parameters[1]
# Base.eltype(branches::Vector{<:Branch}) = typeof(branches).parameters[1].parameters[1]
Base.eltype(tree::Tree) = eltype(tree.branches[1])