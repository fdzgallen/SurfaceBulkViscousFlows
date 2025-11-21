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
τᵈkₒ = 10.0
n  = 40
Δt = 0.001
T  = 0.40
output_frequency = 1

#Rac Coefficients
dᵃ = 200.0  #deactivation
α₀ = 3.00  #basal activation 
Drac = 0.0005 #diffusion rate
αopto = 3 # opto input for Rac 
wrac = π/2.5
χ  = 100.0
χ₀ = -3.0 
rac_total = 3.0  # total amount of Rac

#Rho Coefficients
dᵇ = 200.0 #deactivation
β₀ = 3.00 #basal activation 
Drho = 0.0005 #diffusion rate 
βopto = 3 # opto input for Rho 
σₐ⁰ = 0.4
rho_total = 5.0  # total amount of Rho



# name="SurfaceViscousFlows/velocity_on/rho-rac 5 better friction T=$T dt=$Δt alpha=$α₀ beta=$β₀ da=$dᵃ db= $dᵇ Drac=Drho=$Drac/sig_a=$σₐ⁰ x0=$χ₀ x=$χ sigmarho=1 sigmaR=1 sa=1/"
name="SurfaceViscousFlows/conserved/test lowrac T=$T dt=$Δt alpha=$α₀ beta=$β₀ da=$dᵃ db= $dᵇ Drac=Drho=$Drac/sig_a=$σₐ⁰ rac_t=$rac_total rho_t=$rho_total  x0=$χ₀ x=$χ sigmarho=1 sigmaR=1 sa=1/"
mkpath(name)

GridapPETSc.with() do

  # surface_viscous_flows_axisymmetric(
  #   dᵃ,α₀,Drac,αopto,dᵇ,β₀,Drho,βopto,wrac,
  #   domain,ls,Pe,n,Δt,T,σₐ⁰,χ₀,χ,output_frequency=output_frequency,
  #   writesol=true,initial_density=verification,γᶜ=1.0,τᵈkₒ=τᵈkₒ,
  #   name=name)

  surface_viscous_flows_axisymmetric_conserved(
    dᵃ,α₀,Drac,αopto,dᵇ,β₀,Drho,βopto,wrac,
    domain,ls,Pe,n,Δt,T,rac_total,rho_total,σₐ⁰,χ₀,χ,output_frequency=output_frequency,
    writesol=true,initial_density=verification,γᶜ=1.0,τᵈkₒ=τᵈkₒ,
    name=name)
end