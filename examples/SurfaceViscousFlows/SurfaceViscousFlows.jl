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
T  = 0.40
output_frequency = 1

#Rac Coefficients
dᵃ = 2.0  #deactivation
α₀ = 3.00  #basal activation 
Drac = 0.00005 #diffusion rate
αopto = 20 # opto input for Rac 
wrac=π/2.5

#Rho Coefficients
dᵇ = 2.0 #deactivation
β₀ = 3.00 #basal activation 
Drho = 0.00005 #diffusion rate 
βopto = 0 # opto input for Rho 

name="SurfaceViscousFlows/test/T=$T dt=$Δt alpha=$α₀ beta=$β₀ da=$dᵃ db= $dᵇ Drac=Drho=$Drac/"
mkpath(name)

GridapPETSc.with() do

  surface_viscous_flows_axisymmetric(
    dᵃ,α₀,Drac,αopto,dᵇ,β₀,Drho,βopto,wrac,
    domain,ls,Pe,n,Δt,T,output_frequency=output_frequency,
    writesol=true,initial_density=verification,γᶜ=1.0,τᵈkₒ=τᵈkₒ,
    name=name)

end