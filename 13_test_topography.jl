using Printf
using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBoundary
using Oceananigans.ImmersedBoundaries: solid_node
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Plots
using JLD2
using Oceananigans.ImmersedBoundaries: mask_immersed_field!

function show_mask(grid)

    print("grid = ", grid, "\n")
    c = CenterField(CPU(), grid)
    c .= 1

    mask_immersed_field!(c)

    x, y, z = nodes(c)

    return x, y, z, c
end

grid = RegularRectilinearGrid(size=(16, 8),
                              y=(-1, 1),                        
                              z=(-1, 0),                           
                              topology=(Flat, Periodic, Bounded))

# Gaussian seamount
h0, L = 0.5, 0.25                                        
seamount(x, y, z) = z < - 1 + h0*exp(-y^2/L^2)  
grid_with_seamount = ImmersedBoundaryGrid(grid, GridFittedBoundary(seamount))

###

Ny = grid_with_seamount.grid.Ny
topography_index = zeros(Int64, Ny)
for (iyC, yC) in enumerate(grid_with_seamount.grid.yC[1:Ny])
    print("y = ", yC, " z = ", grid_with_seamount.grid.zC[1], " ")

    izC = 1
    test = true
    while test && izC <= grid_with_seamount.grid.Nz
    	  print(izC, " ")
	  test =  solid_node(Center(), Center(), Center(), 1, iyC, izC, grid_with_seamount)
	  izC += 1
    end
    print(" ")
    print(solid_node(Center(), Center(), Center(), 1,iyC,1, grid_with_seamount)," ")
    print(solid_node(Center(), Center(), Center(), 1,iyC,2, grid_with_seamount),"\n")
    topography_index[iyC] = izC - 1 
end

print(topography_index, "\n")

# suppose we have a tracer and we are calling it Theta

Theta_topography_slice = zeros(Ny)
for (iyC, yC) in enumerate(grid_with_seamount.grid.yC[1:Ny])
    Theta_topography_slice[iyC] = Theta[1, iyC, topography_index[iyC]]
end

plt = plot(yC[1:Ny], Theta_topography_slice)
display(plt)
