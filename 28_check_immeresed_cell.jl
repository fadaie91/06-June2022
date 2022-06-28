using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom, PartialCellBottom


arch = CPU()
grid =  RectilinearGrid(arch,
                                  size=(160, 20), halo=(3, 3), 
                                  y = (-4, 4),
                                  z = (-1, 0),
                                  topology=(Flat, Periodic, Bounded))


h₀ = 0.1 # bump height
L = 1 # bump width
@inline h(y) = h₀ * exp(- y^2 / L^2)
@inline seamount(x, y) = - 1 + h(y)
                                  
seamount_field = Field{Center, Center, Nothing}(grid)
grid_with_seamount = ImmersedBoundaryGrid(grid, GridFittedBottom(seamount_field.data))
immersed_cell(1, 1, 1, grid_with_seamount.underlying_grid)
