using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC

using SurfaceBulkViscousFlows

domain = (-1.2,1.2,0.0,1.2)

ls = AlgoimCallLevelSetFunction(
  x -> x[1]*x[1] + x[2]*x[2] - 1.0,
  x -> VectorValue( 2.0 * x[1], 2.0 * x[2] ) )

Pe = 30.0
τᵈkₒ = 10.0
n  = 30
Δt = 0.0001
T  = 3.0
output_frequency = 1

GridapPETSc.with() do

  surface_viscous_flows_axisymmetric(
    domain,ls,Pe,n,Δt,T,output_frequency=output_frequency,
    writesol=true,initial_density=verification,γᶜ=1.0,τᵈkₒ=τᵈkₒ,
    name="examples/SurfaceViscousFlows/plt")

end