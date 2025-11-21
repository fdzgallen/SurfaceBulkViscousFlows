using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC

using SurfaceBulkViscousFlows

# CAN USE IN TERMINAL BEFORE RUNNING THIS TO ACCELERATE RUNNING TIME BY MULTITHREADING
# export OMP_NUM_THREADS=1 && export JULIA_NUM_THREADS=8

domain = (-1.2,1.2,0.0,1.2)

ls = AlgoimCallLevelSetFunction(
  x -> x[1]*x[1] + x[2]*x[2] - 1.0,
  x -> VectorValue( 2.0 * x[1], 2.0 * x[2] ) )

Pe = 30.0
τᵈkₒ = 0.1
n  = 30
Δt = 0.0002
T  = 0.20
output_frequency = 1

#Rac Coefficients
dᵃ = 200.0  #deactivation
α₀ = 3.00  #basal activation 
Drac = 0.00005 #diffusion rate
αopto = 5 # opto input for Rac 
wrac = π/2.5
χ  = 100.0
χ₀ = -3.0 

#Rho Coefficients
dᵇ = 200.0 #deactivation
β₀ = 3.00 #basal activation 
Drho = 0.00005 #diffusion rate 
βopto = 0 # opto input for Rho 
σₐ⁰ = 0.5

De=0.00000001

name="SurfaceViscousFlows/turnover/rho-rac 5 better friction T=$T dt=$Δt alpha=$α₀ beta=$β₀ da=$dᵃ db= $dᵇ Drac=Drho=$Drac/sig_a=$σₐ⁰ tk=$τᵈkₒ  x0=$χ₀ x=$χ sigmarho=1 sigmaR=1 sa=1/"
mkpath(name)

GridapPETSc.with() do

  surface_viscous_flows_axisymmetric_turnover(
    dᵃ,α₀,Drac,αopto,dᵇ,β₀,Drho,βopto,wrac, De,
    domain,ls,Pe,n,Δt,T,σₐ⁰,χ₀,χ,output_frequency=output_frequency,
    writesol=true,initial_density=verification,γᶜ=1.0,τᵈkₒ=τᵈkₒ,
    name=name)

end