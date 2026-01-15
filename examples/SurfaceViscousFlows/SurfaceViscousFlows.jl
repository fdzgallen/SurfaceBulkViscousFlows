using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC

using SurfaceBulkViscousFlows

# CAN USE IN TERMINAL BEFORE RUNNING THIS TO ACCELERATE RUNNING TIME BY MULTITHREADING
# export OMP_NUM_THREADS=1 && export JULIA_NUM_THREADS=8

domain = (-1.2,0.9,0.0,1.2)
#domain = (-12.0,8.0,0.0,12.0)

R = 1.0 # [um] Cell size length

ls = AlgoimCallLevelSetFunction(
  x -> x[1]*x[1] + x[2]*x[2] - R*R,
  x -> VectorValue( 2.0 * x[1], 2.0 * x[2] ) )

Pe = 30.0
τᵈkₒ = 10.0
n  = 40
Δt = 0.001
T  = 0.900
output_frequency = 1

#Rac Coefficients
dᵃ = 40.0  #deactivation
α₀ = 1.00  #basal activation 
α = 2.00  #  activation 
Drac = 0.01 #diffusion rate
αopto = 1.5 #15 # opto input for Rac 
wrac = π/2 # width of the opto activation
χ  =  0.36 #1.0 #100.0
χ₀ = 0.0 #-3.0 
rac_total = 3.0  # total amount of Rac

#Rho Coefficients
dᵇ = 40.0 #deactivation
β₀ = 1.00 #basal activation 
β = 2.00 #  activation 
Drho = 0.01  #diffusion rate 
βopto = 0#5 # opto input for Rho 
rho_total = 5.0  # total amount of Rho

koff = 1.4*1000  # ezrin koff [1/s] toff=0.7s Fritzsche et al
kon  = 5.0*1000  # ezrin kon [1/s] ton =0.2s Fritzsche et al
D =  0.3/0.036 #30.0   # diffusion of the membrane [um^2/s] 0.003  Fritzsche et al
M0 = 10.0 #total amount of ezrin
#in case one wants nonlinear deactication if the system would be unstable
λᵇ = 0.0011   # for bound ezrin 
λʳᴬ= 0   #for Rho
η = 1000.0  # [pN s/ um] 2D viscosity of the cortex
#λ = sqrt( η*R*π /(χ*R*R*π*π*M0) ) #\bar λ in the suppl material. Adimensional lenthscale
σₐ⁰ = 0.0#.4
sigmaₐ⁰ = 1.0 #basal active tension for the cortex
sigmaρ⁰ = 10.0*sigmaₐ⁰ #rho dependent active tension for the cortex

k = 10.0            # [pN /um] elastic constant for the membrane
rac0 = 0.01 # ADIMENSINAlizes rac switch for protrusion, changes the slope of the threshold function
vCTE = -1.5/0.006 #-0.3 #Scale for the polimerizatio velocity
ten0 = 0.010 #denominator for protrusion dependence on tension

tenth = 0.1 #threshold for tension to activate rho
sig0 = 0.01  #adimensionalizes the tension term for rho
MCAbth = 2.4 #MCA/ezrin value below which protrusion starts
rho0 = 0.001 #adimensionalizes the MCA concentration term for rac
ξ = 1.0
Λ = sqrt(3)*0.25*k/ξ
M = sqrt(3)*0.25*k/ξ
Λ = trunc(Λ,digits=2)
M = trunc(M,digits=2)

# name="SurfaceViscousFlows/velocity_on/rho-rac 5 better friction T=$T dt=$Δt alpha=$α₀ beta=$β₀ da=$dᵃ db= $dᵇ Drac=Drho=$Drac/sig_a=$σₐ⁰ x0=$χ₀ x=$χ sigmarho=1 sigmaR=1 sa=1/"
#name="SurfaceViscousFlows/conserved/test lowrac T=$T dt=$Δt alpha=$α₀ beta=$β₀ da=$dᵃ db= $dᵇ Drac=Drho=$Drac/sig_a=$σₐ⁰ rac_t=$rac_total rho_t=$rho_total  x0=$χ₀ x=$χ sigmarho=1 sigmaR=1 sa=1/"
name="SurfaceViscousFlows/mechanochemical/T=$T dt=$Δt alpha=$α₀ beta=$β₀ da=$dᵃ db= $dᵇ Drac=Drho=$Drac/dimmy vC=$vCTE k=$k sig_a=$sigmaρ⁰ a=$α b=$β  x0=$χ₀ x=$χ/"
mkpath(name)

GridapPETSc.with() do

  # surface_viscous_flows_axisymmetric(
  #   dᵃ,α₀,Drac,αopto,dᵇ,β₀,Drho,βopto,wrac,
  #   domain,ls,Pe,n,Δt,T,σₐ⁰,χ₀,χ,output_frequency=output_frequency,
  #   writesol=true,initial_density=verification,γᶜ=1.0,τᵈkₒ=τᵈkₒ,
  #   name=name)

  # surface_viscous_flows_axisymmetric_conserved(
  #   dᵃ,α₀,Drac,αopto,dᵇ,β₀,Drho,βopto,wrac,
  #   domain,ls,Pe,n,Δt,T,rac_total,rho_total,σₐ⁰,χ₀,χ,output_frequency=output_frequency,
  #   writesol=true,initial_density=verification,γᶜ=1.0,τᵈkₒ=τᵈkₒ,
  #   name=name)

  run_mechanochemical_axisymmetric_vector(koff,kon,M0,D ,λᵇ,λʳᴬ,rac0,rho0,ten0,vCTE,tenth,sig0,MCAbth,
  dᵃ,α₀,Drac,αopto,dᵇ,β₀,Drho,βopto,wrac,sigmaₐ⁰, sigmaρ⁰ , Λ,M, k,
  α,β,
  domain,ls,Pe,n,Δt,T,σₐ⁰,χ₀,χ,output_frequency=output_frequency,
  writesol=true,initial_density=verification,γᶜ=1.0,τᵈkₒ=τᵈkₒ,
  name=name)
end