using Gridap
using GridapEmbedded
using GridapPETSc
using GridapPETSc: PETSC

using SurfaceBulkViscousFlows

# CAN USE IN TERMINAL BEFORE RUNNING THIS TO ACCELERATE RUNNING TIME BY MULTITHREADING
# export OMP_NUM_THREADS=1 && export JULIA_NUM_THREADS=8

#domain = (-1.2,0.9,0.0,1.2)
domain = (-1.2*6.0,0.9*6.0,0.0,1.2*6.0)

R = 6.0 # [um] Cell size length

ls = AlgoimCallLevelSetFunction(
  x -> x[1]*x[1] + x[2]*x[2] - R*R,
  x -> VectorValue( 2.0 * x[1], 2.0 * x[2] ) )

  ###############CONSTANT MECHANICAL PARAMETERS ##############

# # 1. Mecánica
# mech = SurfaceBulkViscousFlows.MechanicalParams(
#     η = 10000.0, k = 50.0, χ = 1.0, χ₀ = 0.0,
#     sigmaₐ⁰ = 1.0, sigmaρ⁰ = 200.0, Λ = 4.3, M = 4.3, R = 6.0
# )

# # 2. Cinética Química
# kin = SurfaceBulkViscousFlows.KineticParams(
#     koff = 1.4, kon = 5.0, M0 = 10.0, D = 0.3,
#     dᵃ = 0.04, dᵇ = 0.04, λᵇ = 0.0, λʳᴬ = 0.0,
#     Drac = 0.0, Drho = 0.0, α₀ = 1.0, β₀ = 1.0, wrac = 4.0
# )

# # 3. Acoplamiento y Optogenética
# coup = SurfaceBulkViscousFlows.CouplingParams(
#     rac0 = 0.01, rho0 = 0.01, ten0 = 10.0, vCTE = - 0.15,
#     tenth = 0.10, sig0 = 1.0, MCAbth = 0.41,
#     α = 1.0, β = 1.0, αopto = 1.0, βopto = 0.0
# )

# # 4. Control de Simulación
# control = SurfaceBulkViscousFlows.SimControl(
#     domain = (-1.2*6.0 , 0.9*6.0 , 0.0 , 1.2*6.0),
#     ls = ls, n = 100, Δt = 1.0, T = 1000.0, order = 2,
#     output_frequency = 1, γᶜ = 1.0, τᵈkₒ = 10.0
# )
ξ = 1.00#2
k = 40.0
# 1. Mecánica
mech = SurfaceBulkViscousFlows.MechanicalParams(
    η = 10000.0, χ = 2.0, χ₀ = 0.0, sigmaₐ⁰ = 100.0, sigmaρ⁰ = 2000.0, S = trunc( sqrt(3)*k*(1 - 1/ξ) , digits=2 ),
    Λ = trunc(sqrt(3)*0.25*k/ξ,digits=2), M = trunc(sqrt(3)*0.25*k/ξ,digits=2), R = 6.0
)
@show mech.Λ
@show mech.M
# 2. Cinética Química
kin = SurfaceBulkViscousFlows.KineticParams(
    koff = 1.4, kon = 5.0, M0 = 1.0, D = 0.3,
    dᵃ = 0.04, dᵇ = 0.04, λᵇ = 0.0, λʳᴬ = 0.0,
    Drac = 0.05, Drho = 0.05, α₀ = 1.0, β₀ = 1.0, wrac = 6.0
)

# 3. Acoplamiento y Optogenética
coup = SurfaceBulkViscousFlows.CouplingParams(
    rac0 = 0.01, rho0 = 0.001, ten0 = 0.02, vCTE = - 0.3,
    tenth = 0.025, sig0 = 0.01, MCAbth = 0.035,
    α = 1.0, β = 1.0, αopto = 2.0, βopto = 0.0
)

# 4. Control de Simulación
control = SurfaceBulkViscousFlows.SimControl(
    domain = (-1.2*6.0 , 0.9*6.0 , 0.0 , 1.2*6.0),
    ls = ls, n = 40, Δt = 1.0, T = 350.0, order = 2,
    output_frequency = 1, γᶜ = 1.0, τᵈkₒ = 10.0
)

 
# Calculate lengthscale
len = trunc((sqrt( mech.η*mech.R*pi /(mech.χ*mech.R*mech.R*kin.M0*pi*pi) )),digits=2)

println("Lengthscale = $len")

# name="SurfaceViscousFlows/velocity_on/rho-rac 5 better friction T=$T dt=$Δt alpha=$α₀ beta=$β₀ da=$dᵃ db= $dᵇ Drac=Drho=$Drac/sig_a=$σₐ⁰ x0=$χ₀ x=$χ sigmarho=1 sigmaR=1 sa=1/"
#name="SurfaceViscousFlows/conserved/test lowrac T=$T dt=$Δt alpha=$α₀ beta=$β₀ da=$dᵃ db= $dᵇ Drac=Drho=$Drac/sig_a=$σₐ⁰ rac_t=$rac_total rho_t=$rho_total  x0=$χ₀ x=$χ sigmarho=1 sigmaR=1 sa=1/"
#name="SurfaceViscousFlows/mechanochemical/T=$T dt=$Δt alpha=$α₀ beta=$β₀ da=$dᵃ db= $dᵇ Drac=Drho=$Drac/Dimensional vC=$vCTE k=$k sig_a=$sigmaρ⁰ a=$α b=$β  x0=$χ₀ x=$χ/"
name = "SurfaceViscousFlows/mechanochemical/Final/" * "T=$(control.T)_dt=$(control.Δt)_ n=$(control.n)/D=$(kin.D)_" *
       "eta=$(mech.η)_vC=$(coup.vCTE) S=$(mech.S)_M=$(mech.M)_a0=$(kin.α₀) b0=$(kin.β₀) wrac=$(kin.wrac) Drac=$(kin.Drac)" *
       "/rac0=$(coup.rac0) koff=$(kin.koff) aopto = $(coup.αopto) da=$(kin.dᵃ) db=$(kin.dᵇ)_sig=$(mech.sigmaρ⁰)_x=$(mech.χ)  L=$len a=$(coup.α) b=$(coup.β) sig0=$(coup.sig0)_tenth=$(coup.tenth)_MCAbth=$(coup.MCAbth)_rho0=$(coup.rho0) ten0=$(coup.ten0)/" 
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

  # run_mechanochemical_axisymmetric_vector(koff,kon,M0,D ,λᵇ,λʳᴬ,rac0,rho0,ten0,vCTE,tenth,sig0,MCAbth,
  # dᵃ,α₀,Drac,αopto,dᵇ,β₀,Drho,βopto,wrac,sigmaₐ⁰, sigmaρ⁰ , Λ,M, k, η,
  # α,β,
  # domain,ls,n,Δt,T,χ₀,χ,output_frequency=output_frequency,
  # writesol=true,initial_density=verification,γᶜ=1.0,τᵈkₒ=τᵈkₒ,
  # name=name)
  run_mechanochemical_axisymmetric_vector(mech, kin, coup, control, name=name, writesol=true)
end