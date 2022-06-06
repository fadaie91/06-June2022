using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom, PartialCellBottom
using Printf
using Plots
using Oceananigans.ImmersedBoundaries: mask_immersed_field!

arch = CPU()

function show_mask(grid)

    c = CenterField(grid)
    c .= 1

    mask_immersed_field!(c)

    x, y, z = nodes(c)

    return x, y, z, c
end

underlying_grid = RectilinearGrid(arch,
                                         size=(640, 80), halo=(3, 3),
                                         y = (-4, 4),
                                         z = (-1, 0),
                                         topology=(Flat, Periodic, Bounded))

# A bump
h₀ = 0.1 # bump height
L = 1 # bump width

@inline h(y) = h₀ * exp(- y^2 / L^2)
@inline seamount(x, y) = - 1 + h(y)

seamount_field = Field{Center, Center, Nothing}(underlying_grid)
set!(seamount_field, seamount)
fill_halo_regions!(seamount_field)
grid_with_seamount = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(seamount_field.data))

x, y, z, cc = show_mask(grid_with_seamount)

plt_full = heatmap(y, z, interior(cc)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region_Full_cell")
savefig(plt_full,"mask_fullcell.png")


grid_with_seamount = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(seamount_field.data,minimum_fractional_Δz=1))

x, y, z, cc = show_mask(grid_with_seamount)

plt_partial = heatmap(y, z, interior(cc)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region_Partial_cell")
savefig(plt_partial,"mask_partialcell.png")
