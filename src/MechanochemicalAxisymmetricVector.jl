# Standard code to run a simulation using the mechanochemical model presented in article "INSERT DOI"
# The system will start with flat Rho and Rac at an steady state
# Different mechanical parametes can eb changed in the main section of the code.
# One can work with pure local inhibition by changing the following mechanical parameters:
# coup.vCTE=0   coup.α=0   coup.β=0   σₐ₀=0
# Written by Andreu F Gallen (working in Turlier lab) and Eric Neiva, in collaboration with Orion Weiner's lab
include("Plots_RhoRacA.jl") 

@kwdef struct MechanicalParams
    η::Real           # Viscosidad
    χ::Real           # Fricción
    χ₀::Real          # Fricción basal
    sigmaₐ⁰::Real     # Tensión activa basal
    sigmaρ⁰::Real     # Tensión dependiente de Rho
    S::Real           # Parámetro geométrico/estabilización
    Λ::Real           # Parámetro geométrico/estabilización
    M::Real           # Parámetro geométrico/estabilización
    R::Real           # Radio de la célula
end

@kwdef struct KineticParams
    koff::Real; kon::Real; M0::Real; D::Real
    dᵃ::Real; dᵇ::Real; λᵇ::Real; λʳᴬ::Real
    Drac::Real; Drho::Real; α₀::Real; β₀::Real; wrac::Real
end

@kwdef struct CouplingParams
    rac0::Real; rho0::Real; ten0::Real; vCTE::Real
    tenth::Real; sig0::Real; MCAbth::Real
    α::Real; β::Real; αopto::Real; βopto::Real
end

@kwdef struct SimControl
    domain::Tuple
    ls::Any           # El objeto LevelSet
    n::Int            # Partición
    Δt::Real
    T::Real
    order::Int
    output_frequency::Int
    γᶜ::Real
    τᵈkₒ::Real
end


function plotting(ylab,po,pPNG,i)
  plot(po)
  xlabel!("ξ[μm]")
  ylabel!(ylab)
  savefig(pPNG*"$ylab"*i*".png")
end

function conservation(sMCAu,sMCAb,Minitial)
  return 40.0*(Minitial-(sMCAu+sMCAb))
end
 
function threshold(x,x₀,xth)
  return  (0.5 * (tanh∘(x/x₀ - xth/x₀)+1))
end
function threshold2(x,x₀,xth)
  return  (0.5 * (tanh.(x/x₀ .- xth/x₀).+1))
end

function sinθ(x)
  return x[2]/norm(x)
end

function cotθ(x)
  return  x[1]/x[2]  #  cosθ(x)/sinθ(x)
end

function cotθ2(x)
  return  cotθ(x)*cotθ(x)
end

function cotθ3(x)
  return  cotθ(x)*cotθ(x)*cotθ(x)
end

function cscθ2(x)
  return  1 + cotθ2(x) 
end


function MCA_bound_unbound_weak_forms(Δt,kon,koff,λᵇ,λ,R2,D,nΓ,dΓ)
    # Now for MAC bound
  # mass term for the temporal evolution MCA_b 
  τ = TensorValue(0.0,-1.0, 1.0, 0.0) ⋅ nΓ    # vector tangente

  mMCA(Δt,MCA_b,w) = ∫( ( (MCA_b*w)/Δt )*y )dΓ 
  aMCAb(MCA_b,v,w) = ∫( ( (MCA_b*w)/Δt )*y )dΓ + ∫( ( 0.0005 * (∇ᵈ(MCA_b,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ +  #diffusion term for stabilisation purposes only
    ∫( ( koff * ( MCA_b * w ) )*y )dΓ + 
    ∫( ( w * ( (v ⋅ τ) * (∇ᵈ(MCA_b,nΓ)⋅τ)))*y )dΓ  +  
    ∫( w * ( MCA_b * ( (τ ⋅ ∇ᵈ(v,nΓ)) ⋅ τ  ) ) * y )dΓ  
  bMCAb(w,MCA_u,MCAb_old,λ) = ∫( ( kon * ( MCA_u * w ) + (MCAb_old*w)/Δt )*y )dΓ  +
    ∫( ( λ/(π*R2) * ( kon / (kon+koff) )* w )*y )dΓ

  #Now for MCA unbound
  # mass term for the temporal evolution MCA_b0
  aMCAu(MCA_u,x,x_old,w) = ∫( ( (MCA_u*w)/Δt )*y )dΓ + ∫( ( D * ( ∇ᵈ(MCA_u,nΓ)⋅∇ᵈ(w,nΓ) ) )*y )dΓ + 
    ∫( ( kon * ( MCA_u * w ) )*y )dΓ + 
    ∫( ( w * ( ( (x-x_old)⋅τ / Δt ) * ( ∇ᵈ(MCA_u,nΓ)⋅τ ) + 
    MCA_u * ( ( τ ⋅ ∇ᵈ(x,nΓ)⋅τ  ) - ( τ ⋅ ∇ᵈ(x_old,nΓ)⋅τ ) ) / Δt ) )*y )dΓ  #+ 
    #∫( ( λᵇ * ( ( MCA_u*MCA_u*MCA_u ) * w ) )*y )dΓ
  bMCAu(w,MCA_b,MCAu_old,λ) = ∫( ( koff * (MCA_b*w))*y )dΓ + ∫( ( (MCAu_old*w)/Δt )*y )dΓ  +
    ∫( ( λ/(π*R2) * (koff/(kon+koff))* w )*y )dΓ # mMCA(control.Δt,MCAu_old,w) 

  mMCA,aMCAb,bMCAb,aMCAu,bMCAu
end
 
#function to run a single simulation with a few given parameters
function run_mechanochemical_axisymmetric_vector(
    mech::MechanicalParams,
    kin::KineticParams,
    coup::CouplingParams,
    control::SimControl; 
    name::String, 
    writesol::Bool = true
)

    activity::Function = unit_activity_axisymmetric
    redistance_frequency::Int = 1
  # Background geometry
  cells = (control.n,div(control.n,2))
  h = (control.domain[2]-control.domain[1])/control.n
  bgmodel = CartesianDiscreteModel(control.domain,cells)
  Ω = Triangulation(bgmodel)

  # Buffer of active model and integration objects
  degree = control.order < 3 ? 3 : 2*control.order
  buffer = Ref{Any}(( Ωᶜ  = nothing, dΩᶜ = nothing,
                      dΓ  = nothing, nΓ  = nothing,
                      φ₋  = nothing, cp₋ = nothing, 
                      t   = nothing, Vbg = nothing ))

  function update_buffer!(i,t,dt,v₋₂,mv₋₂)

    if buffer[].t == t
      return true
    else

      Ωᶜ    = buffer[].Ωᶜ
      dΩᶜ   = buffer[].dΩᶜ
      dΓ    = buffer[].dΓ
      nΓ    = buffer[].nΓ
      cp₋   = buffer[].cp₋
      φ₋    = buffer[].φ₋
      t     = buffer[].t
      Vbg   = buffer[].Vbg

      if buffer[].Ωᶜ === nothing
        Vbg = TestFESpace(Ω,ReferenceFE(lagrangian,Float64,control.order))
        _φ₋ = interpolate_everywhere(control.ls.φ,Vbg)
      else
        cp₋₂ = buffer[].cp₋
        φ₋₂  = buffer[].φ₋
        __φ = get_free_dof_values(φ₋₂.φ)
        Ωⱽ  = get_triangulation(Vbg)
        _ϕ₋ = compute_normal_displacement(cp₋₂,φ₋₂,v₋₂,dt,Ωⱽ)
        ϕ₋  = __φ - _ϕ₋
        _φ₋ = FEFunction(Vbg,ϕ₋)
      end

      # Current time level set
      φ₋  = AlgoimCallLevelSetFunction(_φ₋,∇(_φ₋))
      ( i % redistance_frequency == 0 ) && begin
        _φ₋ = compute_distance_fe_function(bgmodel,Vbg,φ₋,control.order,cppdegree=3)
        φ₋  = AlgoimCallLevelSetFunction(_φ₋,∇(_φ₋))
      end

      cp₋ = compute_closest_point_projections(Vbg,φ₋,control.order,
              cppdegree=3,trim=true,limitstol=1.0e-2)

      # Current time surface and bulk measures
      squad = Quadrature(algoim,φ₋,degree,phase=CUT)
      s_cell_quad,is_c₋ = CellQuadratureAndActiveMask(bgmodel,squad)

      # Surface narrow-band triangulation
      δ₋ = 2.0 * mv₋₂ * dt
      _,is_nᶜ = narrow_band_triangulation(Ω,_φ₋,Vbg,is_c₋,δ₋)

      Ωᶜ,dΓ = TriangulationAndMeasure(Ω,s_cell_quad,is_nᶜ,is_c₋)

      dΩᶜ = Measure(Ωᶜ,2*control.order)
      nΓ = normal(φ₋,Ω)

      # Update buffer
      buffer[] = ( Ωᶜ=Ωᶜ,dΩᶜ=dΩᶜ,dΓ=dΓ,nΓ=nΓ,cp₋=cp₋,φ₋=φ₋,t=t,Vbg=Vbg )
      return true

    end

  end

  # Reference FEs
  N = num_dims(bgmodel)
  reffeʷ = ReferenceFE(lagrangian,VectorValue{N,Float64},control.order-1)
  reffeᵉ = ReferenceFE(lagrangian,Float64,control.order-1)

  function update_all!(i::Int,t::Real,dt::Real,disp,val::Real)

    #update_buffer!(i,t,dt,disp,val)

    # Triangulations and aggregates
    Ωᶜ = buffer[].Ωᶜ

    # Measures and normal
    dΩᶜ = buffer[].dΩᶜ
    dΓ  = buffer[].dΓ
    nΓ  = buffer[].nΓ
    φ   = buffer[].φ₋

    # Test FE spaces

    τ = TensorValue(0.0,-1.0, 1.0, 0.0) ⋅ nΓ    # vector tangente 
    # u-surface
    Vʷ = TestFESpace(Ωᶜ,reffeʷ,dirichlet_tags=[5,8]) # vector space test FE space
    UXʷ = TrialFESpace(Vʷ,[p->VectorValue(0,x₀),p->-1.0*xₗ*τ(p)]) # vector space trial FE space for membrane displacement
    UVʷ = TrialFESpace(Vʷ,[p->VectorValue(0,v₀),p->-1.0*vₗ*τ(p)]) # vector space trial FE space for cortex velocity
    
    # e-surface
    Vᵉ = TestFESpace(Ωᶜ,reffeᵉ)
    # Rac-Rho surface
    Vᴿ = TestFESpace(Ωᶜ,reffeᵉ)
    # Lagrange multipliers
    Vˡ = ConstantFESpace(bgmodel)

    # Trial FE spaces
    Uʷ = TrialFESpace(Vʷ)
    Uᵉ = TrialFESpace(Vᵉ)
    Uᴿ = TrialFESpace(Vᴿ)
    Uˡ = TrialFESpace(Vˡ)

    # Multifield FE spaces
    Yᵛ = Vʷ
    Xᵛ = MultiFieldFESpace([Uʷ,Uˡ])
    UXᵛ = MultiFieldFESpace([UXʷ,Uˡ])
    UVᵛ = UVʷ

    # Space to create homogeneous perturbation  
    # of constant concentration myosin field
    # TO-DO: To be deleted if needed
    Yʳ = MultiFieldFESpace([Vᵉ,Vˡ])
    Xʳ = MultiFieldFESpace([Uᵉ,Uˡ])



    UXʷ,UVʷ,Vʷ, UXᵛ,UVᵛ,Xᵛ,Yᵛ,Xʳ,Yʳ,Uᵉ,Vᵉ,Vᴿ,Uᴿ,Ωᶜ,dΩᶜ,dΓ,nΓ,φ

  end

 # Lets make output folders
  pVTU="./output/"*name*"VTU/"
  mkpath(pVTU)
  pPNG="./output/"*name
  mkpath(pPNG)
  mkpath(pPNG*"Rac_time/") 
  mkpath(pPNG*"Rho_time/") 
  mkpath(pPNG*"MCAb_time/") 
  mkpath(pPNG*"Rac_time_initial/") 
  mkpath(pPNG*"Rho_time_initial/") 
 
  # Lets copy the code in the output folder to be able to check code used for each simulation
  cp(@__FILE__, pPNG*split(@__FILE__, "/")[end],force=true)
  cp("./src/WeakForms.jl", pPNG*"WeakForms.jl",force=true)
  cp("./examples/SurfaceViscousFlows/SurfaceViscousFlows.jl", pPNG*"SurfaceViscousFlows.jl",force=true)

  # Time discretisation parameters
  t₀ = 0.0 
  u₀ = VectorValue(0.0,0.0)
  m₀ = 2.0
  R2 = mech.R
  nΔt = trunc(Int,control.T/control.Δt+0.5)+1 

  #Starting Boundary conditions
  x₀ = 0.0
  xₗ = 0.0 
  v₀ = 0.0
  vₗ = 0.0

  update_buffer!(0,t₀,control.Δt,u₀,m₀)
  UXʷ,UVʷ,Vʷ,UXᵛ,UVᵛ,Xᵛ,Yᵛ,Xʳ,Yʳ,Uᵉ,Vᵉ,Vᴿ,Uᴿ,Ωᶜ,dΩᶜ,dΓ,nΓ,φ = update_all!(0,t₀,control.Δt,u₀,m₀)

  τ = TensorValue(0.0,-1.0, 1.0, 0.0) ⋅ nΓ    # vector tangente 

  # *** WEAK FORM PARAMETERS ***
  ξ(e) = 2.0 * e*e / ( 1.0 + e*e )
  # ** u-stabilisation **
  γʷ = control.γᶜ/h
  # ** e-stabilisation **
  γᵉ = control.γᶜ/h


  # Compute initial condition for surface velocity
  _υₕ(x) = VectorValue(0.0,0.0)
  υₕ  = interpolate_everywhere(_υₕ,UVᵛ)
  # Compute initial condition for membrane deformation
  _xₕ(x) = VectorValue(0.0,0.0)
  xₕ  = interpolate_everywhere(_xₕ,UXᵛ[1])
  xₕ_old = xₕ

  Tm = SparseMatrixCSR{0,PetscScalar,PetscInt}
  Tv = Vector{PetscScalar}
  ps = PETScLinearSolver(mykspsetup)

  i = 0
  t = t₀
  
  tol = 1e-8


  #Building the vectors used for introducing opto influence as an increase in coup.α and coup.β
  arclength(x) = R2 * atan(x[2],-x[1]) # Arc length for sphere
  α₀opto(x) = kin.α₀
  β₀opto(x) = kin.β₀

  α₀opto2(x) = kin.α₀ + coup.αopto * exp( -0.5 * ( arclength(x)-π*R2 )^2 / ((kin.wrac)^2) )
  β₀opto2(x) = kin.β₀ + coup.βopto * exp( -0.5 * ( arclength(x))^2       / ((kin.wrac)^2) ) 

  γ₀ = 0.1  / h # TODO: Eric reviews the scaling with h
  γ₀R =  0.1  / h  
  γ₀M =  0.1  / h  
  m₀opto(u,v) = ∫( u*v )dΓ 
  s₀opto(u,v) = ∫( γ₀*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ
  s₀R(u,v) = ∫( γ₀R*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ
  s₀MCA(u,v) = ∫( γ₀M*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ 
  s₀x(υ,μ) = ∫( γʷ * ((nΓ⋅ε(υ))⊙(nΓ⋅ε(μ))) )dΩᶜ
  λ=0.0

  A₀opto(u,v) = m₀opto(u,v) + s₀opto(u,v) 
  bα₀opto(v) = m₀opto(α₀opto,v)
  bβ₀opto(v) = m₀opto(β₀opto,v)
  bα₀opto2(v) = m₀opto(α₀opto2,v)
  bβ₀opto2(v) = m₀opto(β₀opto2,v)

  op_α₀ = AffineFEOperator(A₀opto,bα₀opto,Uᴿ,Vᴿ)
  op_β₀ = AffineFEOperator(A₀opto,bβ₀opto,Uᴿ,Vᴿ)

  α₀v = solve(op_α₀)
  β₀v = solve(op_β₀)
  
  #Rac and Rho initialization.3
  Rₕ = interpolate_everywhere(0.0,Uᴿ) 
  ρₕ = interpolate_everywhere(1.0,Uᴿ) 
  Rₕ_old = Rₕ
  ρₕ_old = ρₕ
  a_R, b_R, a_ρ, b_ρ = rac_rho_weak_forms2(control.Δt,200*kin.dᵃ,200*(kin.dᵇ),kin.Drac,kin.Drho,nΓ,dΓ,coup.α,coup.β)
 
  #ezrin initialization
  mMCA,aMCAb,bMCAb,aMCAu,bMCAu = MCA_bound_unbound_weak_forms(
    control.Δt,kin.kon,kin.koff,kin.λᵇ,λ,R2,kin.D,nΓ,dΓ)

  #we write down the weak form of the membrane equation
  # to do backward eurler for the time evolwe give gridap a so-called Mass Term for a 
  # Weak form of the membrane equation
  #
  # Nonlinear strain rate: 
  # εᴾ(u) = ε(u) + 0.5 * ∇(u)ᵗ⋅∇(u) = εᴾ(u) = ε(u) + εᴺ(u)
  #
  # Terms of the bilinear and linear forms
  # TERM 1. ∫( 2M⋅εᴾ(u):ε(v) )dΓ = ∫( 2M⋅ε(u):ε(v) + 
  #                                   2M⋅εᴺ(u):ε(v) )dΓ 
  #TODO ENSURE IF ITS 2*M (previously) OR M
  aᴹ(M,R,x_old,x,w) = 
    ∫( ( 2* M * ( x⋅w/2 + 
                R*R * ( ∇ᶜ(x,nΓ) ⊙ ∇ᶜ(w,nΓ) ) + 
                (cotθ2) * (x⋅w) ) ) * sinθ )dΓ +
    ∫( ( 2*M/R/2 * ( 
      ( x_old ⋅ x + R*R * ( ∇ᶜ(x_old,nΓ)⊙∇ᶜ(x,nΓ) ) ) * ( R * ∇ᶜ(w,nΓ)⋅τ ) + 
      ( cotθ3 * x_old ) * (x ⋅ w) ) )⋅τ * sinθ )dΓ 
  #
  # TERM 2. ∫( L⋅(tr(εᴾ(u))Id):ε(v) )dΓ = ∫( L⋅tr(ε(u)):ε(v) +  
  #                                          L⋅tr(εᴺ(u)):ε(v) )dΓ
  # Homework: Implement TERM 2
  aᴸ(L,R,x_old,x,w) = 
    ∫( L*( R*R*(∇ᶜ(x,nΓ) ⊙ ∇ᶜ(w,nΓ)) + (cotθ*(x ⋅ ∇ᶜ(w,nΓ)) + cotθ*(w⋅∇ᶜ(x,nΓ)))⋅τ + 
                (cotθ2) * (x⋅w)  ) * sinθ )dΓ +
    ∫( 0.5/R*L*( ( (cscθ2)*(x_old⋅x) + R*R * ∇ᶜ(x,nΓ) ⊙ ∇ᶜ(x_old,nΓ) )*
    ( R*∇ᶜ(w,nΓ)⋅τ + cotθ * w )⋅τ ) * sinθ )dΓ

  bₓ(MCA_b,v,w) = ∫( ( mech.χ*( MCA_b)*v⋅w ) * (mech.R*mech.R) * sinθ )dΓ
  # Preserve mass term for Backward Euler time integration
  m(MCA_b,Δt,x,w) = ∫( ( (mech.χ*MCA_b) * (x⋅w) / Δt ) * (mech.R*mech.R) * sinθ )dΓ
 

  # ** weak tangentiality **
  bo = 10.0 / ((2/40)^2)
  wt(x,w) = ∫( bo*((x⋅nΓ)*(w⋅nΓ)) )dΓ

  mυ(Δt,v,w) = ∫( ( (v⋅w) / Δt )*y )dΓ
  #SOLVING X AT t=0 
  # A(x,w) = m(0,control.Δt,x,w) + a(mech.Λ,mech.M,R2,xₕ,x,w) + s₀x(x,w)
  # B(w) = m(0,control.Δt,xₕ,w) + bₓ(0,υₕ,w)   
  # op_x = AffineFEOperator(A,B,UXᵛ,Yᵛ)
  # xₕ = solve(op_x) 


  #MCA initialization
  uh_MCAb = interpolate_everywhere(kin.kon*kin.M0/(π*R2)/(kin.koff+kin.kon),Uᴿ)
  uh_MCAb_old = uh_MCAb
  uh_MCAu = interpolate_everywhere(kin.koff*kin.M0/(π*R2)/(kin.koff+kin.kon),Uᴿ)
  uh_MCAu_old = uh_MCAu
  sum_uh_MCAu = ∑(∫(uh_MCAu)dΓ)
  sum_uh_MCAb = ∑(∫(uh_MCAb)dΓ)
  Minitial = sum_uh_MCAu + sum_uh_MCAb 
  λ = conservation(sum_uh_MCAu,sum_uh_MCAb,Minitial)

  I = TensorValue(1.0,0.0,0.0,1.0)

 # Computing tension  
  mten(u,v) = ∫( (u*v)*y )dΓ
  N(u) = mech.S * I + mech.Λ * tr(εᶜ(u,nΓ)) * I + mech.M * (εᶜ(u,nΓ))
  # @show N
  # tNt = τ ⋅ (N ⋅ τ)
  # @show tNt
  ∂u(u) = mech.R*(τ ⋅ (∇ᶜ(u,nΓ))⋅ τ)
  mten2(u,v) = ∫(( ( τ ⋅ (N(u) ⋅ τ)  )*v)*y )dΓ # ∫(( ( mech.S + (mech.Λ+mech.M)*(( ∂u(u) + (u ⋅ τ) * cotθ)/mech.R + ( ∂u(u)*∂u(u) + (u⋅u)*cscθ2)/(mech.R*mech.R))   )*v)*y )dΓ
  #mten2(u,v) = ∫(( ( mech.S + (mech.Λ+mech.M)*( ∂u(u) - (u ⋅ τ) * cotθ)/mech.R   )*v)*y )dΓ
  sten(u,v) = ∫( 10*γ₀*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ

  Aten(u,v) = mten(u,v) + sten(u,v) 
  bten(v) = mten2(xₕ,v)
  op_ten = AffineFEOperator(Aten,bten,Uᴿ,Vᴿ)
  
  #_trN = ( 1 + (mech.Λ+mech.M)*( mech.R*(τ ⋅ (∇ᵈ(xₕ,nΓ))⋅ τ) + (xₕ ⋅ τ) * cotθ)/mech.R + ( mech.R*(τ ⋅ (∇ᵈ(xₕ,nΓ))⋅ τ)*mech.R*(τ ⋅ (∇ᵈ(xₕ,nΓ))⋅ τ) + (xₕ⋅xₕ)*cscθ2)/(mech.R*mech.R) )
  ten = solve(op_ten)

  #SOLVE Rac AT t=0
  Arac(rac,w) = a_R(rac,w,υₕ) + s₀R(rac,w)
  Brac(w) = b_R(w,ρₕ,α₀v,Rₕ_old,uh_MCAb ,coup.MCAbth,coup.rho0) #(w,rho,α₀v,rac_old)
  op_rac = AffineFEOperator(Arac,Brac,Uᴿ,Vᴿ)
  Rₕ = solve(op_rac)
  Rₕ_old = Rₕ
  #SOLVE Rho AT t=0
  Arho(rho,w) = a_ρ(rho,w,υₕ) + s₀R(rho,w)
  Brho(w) = b_ρ(w,Rₕ,β₀v,ρₕ_old,ten,coup.sig0,coup.tenth) #(w,rac,β₀v,rho_old)
  op_rho = AffineFEOperator(Arho,Brho,Uᴿ,Vᴿ)
  ρₕ = solve(op_rho) 
  ρₕ_old = ρₕ
  
  #SOLVE MCA_b AT t=0 vien initial velocity zero
  AMCAb(MCA_b,w) = aMCAb(MCA_b,υₕ,w) + s₀MCA(MCA_b,w)
  BMCAb(w) = bMCAb(w,uh_MCAu,uh_MCAb_old,λ) 

  #SOLVE MCA_u AT t=0
  AMCAu(MCA_u,w) = aMCAu(MCA_u,xₕ,xₕ_old,w) + s₀MCA(MCA_u,w)
  BMCAu(w) = bMCAu(w,uh_MCAb,uh_MCAu_old,λ) 
  
  # Extract quadrature points and arc length array at every cell
  xΓ = dΓ.quad.cell_point.values
  xΓ = lazy_map(Reindex(xΓ),dΓ.quad.cell_point.ptrs)     # 2D array of xΓ (1 array per cell)
  alenΓ = lazy_map(Broadcasting(x->atan(x[2],-x[1])),xΓ) # Following cell order
  flat_xΓ = vcat(xΓ...)                                  # 1D "flattened" array of xΓ
  flat_alenΓ = vcat(alenΓ...)
  num_qpoints = length(flat_xΓ)
  perm=sortperm(flat_alenΓ) # Permutation to order by increasing arclength
  flat_alenΓ = R2*flat_alenΓ[perm]

  
  ract = zeros(nΔt,num_qpoints)
  rhot = zeros(nΔt,num_qpoints)
  MCAbt = zeros(nΔt,num_qpoints)
  σₐt = zeros(nΔt,num_qpoints)
  χt = zeros(nΔt,num_qpoints)
  vt = zeros(nΔt,num_qpoints)
  vnt = zeros(nΔt,num_qpoints)
  xt = zeros(nΔt,num_qpoints)
  tent = zeros(nΔt,num_qpoints)
  _vt = vcat(lazy_map(υₕ ⋅ τ,xΓ)...) 
  vt[1,:] = _vt[perm] 
  _vnt = vcat(lazy_map(υₕ ⋅ nΓ,xΓ)...) 
  vnt[1,:] = _vnt[perm] 
  _xt = vcat(lazy_map(xₕ ⋅ τ,xΓ)...) 
  xt[1,:] = _xt[perm] 
  MCAbt[1,:] = vcat(lazy_map(uh_MCAb,xΓ)...)[perm]
  tent[1,:] = vcat(lazy_map(ten,xΓ)...)[perm] 


  function sigmaₐ(ρ,R)
      sigmaₐ = mech.sigmaₐ⁰ .+ mech.sigmaρ⁰ * ρ #.- sigmaR⁰ * R
      #sigmaₐ > 0 ? sigmaₐ : zero(typeof(sigmaₐ))
  end

  χR(R) = (mech.χ₀.+mech.χ*uh_MCAb)
 

  for ti in 1:100 
    op_rho = AffineFEOperator(Arho,Brho,Uᴿ,Vᴿ)
    ρₕ = solve(op_rho)
    ρₕ_old = ρₕ
    op_rac = AffineFEOperator(Arac,Brac,Uᴿ,Vᴿ)
    Rₕ = solve(op_rac)
    Rₕ_old = Rₕ
    sum_R = ∑(∫(Rₕ)dΓ)
    sum_ρ = ∑(∫(ρₕ)dΓ)
    ractt = vcat(lazy_map(Rₕ,xΓ)...)
    rhott = vcat(lazy_map(ρₕ,xΓ)...) 
    ractt[:] = ractt[perm]
    rhott[:] = rhott[perm]
    plotting("rac",ractt[:],pPNG*"Rac_time_initial/","$ti")
    plotting("rho",rhott[:],pPNG*"Rho_time_initial/","$ti")
  end

  msₕ = get_maximum_magnitude_with_dirichlet(υₕ)
  a_R, b_R, a_ρ, b_ρ = rac_rho_weak_forms2(control.Δt,kin.dᵃ,kin.dᵇ,kin.Drac,kin.Drho,nΓ,dΓ,coup.α,coup.β)
  Arac(rac,w) = a_R(rac,w,υₕ) + s₀R(rac,w)

  while t < control.T + tol
    if i > 50
      op_α₀ = AffineFEOperator(A₀opto,bα₀opto2,Uᴿ,Vᴿ)
      op_β₀ = AffineFEOperator(A₀opto,bβ₀opto2,Uᴿ,Vᴿ)
      α₀v = solve(op_α₀)
      β₀v = solve(op_β₀)
    end
    if i > 150
      op_α₀ = AffineFEOperator(A₀opto,bα₀opto,Uᴿ,Vᴿ)
      op_β₀ = AffineFEOperator(A₀opto,bβ₀opto,Uᴿ,Vᴿ)
      α₀v = solve(op_α₀)
      β₀v = solve(op_β₀)
    end
    
    i1 = ∑(∫(uh_MCAb)dΓ)
    i2 = ∑(∫(uh_MCAu)dΓ)
    λ = conservation(i2,i1,Minitial)

    @info "Time step $i, time $(trunc(t, digits=4)) and time step $(control.Δt)"
 
    
    aᵛ,bᵛ = cortical_flow_problem_mechanochemical_axisymmetric_dimensional(xₕ,xₕ_old, control.Δt,mech.η,
         ρₕ,uh_MCAb,dΩᶜ,dΓ,nΓ,γʷ,mech.χ,mech.χ₀,activity,mech.sigmaₐ⁰,  mech.sigmaρ⁰)
    op = AffineFEOperator(aᵛ,bᵛ,UVᵛ,Yᵛ)
    υₕ = solve(op)
    υₕtan = to_tangent_vector(υₕ,nΓ) #υₕ⋅(TensorValue(0.0,-1.0,1.0,0.0)⋅nΓ)#⋅(VectorValue(0.0,-1.0,1.0,0.0)⋅nΓ)

    #msₕ = get_maximum_magnitude_with_dirichlet(υₕ)

    aˣ(x,w) = m(uh_MCAb,control.Δt,x,w) + aᴸ(mech.Λ,R2,xₕ,x,w) + aᴹ(mech.M,R2,xₕ,x,w) + s₀x(x,w) + wt(x,w) #a(mech.Λ,mech.M,R2,xₕ,x,w) 
    bˣ(w) = m(uh_MCAb,control.Δt,xₕ,w) + bₓ(uh_MCAb,υₕ,w)    

    op_x= AffineFEOperator(aˣ,bˣ,UXʷ,Vʷ)
    xₕ = solve(op_x) 
    op_ten = AffineFEOperator(Aten,bten,Uᴿ,Vᴿ) #Utʷ,Vᵗ) 
    # N = mech.S * I + mech.Λ * tr(εᶜ(xₕ,nΓ)) * I + mech.M * (εᶜ(xₕ,nΓ))
    # tNt = τ ⋅ (N ⋅ τ)
    ten = solve(op_ten)
   # ten = ∫(( mech.k*( τ ⋅ ∇ᵈ(xₕ,nΓ) ⋅ τ ))*y )dΓ
    xₕ_old =  xₕ
 
    i = i + 1
    t = t + control.Δt

    writesol && postprocess_all_with_tangent(φ,dΩᶜ.quad.trian,
      Rₕ,ρₕ,xₕ,υₕ,υₕtan,uh_MCAb,ten,i=i,of=control.output_frequency,name=pVTU)

    tent[i,:] = vcat(lazy_map(ten,xΓ)...)[perm]
    #boundary conSSditions for v and x 
    Rₕaux = vcat(lazy_map(Rₕ,xΓ)...)[perm]   
    _ten = vcat(lazy_map(ten,xΓ)...)[perm]
    _ten = _ten[end]*threshold2( _ten[end] , 0.0001 , 0.0 )
    vₗ = coup.vCTE * threshold2( Rₕaux[end] , coup.rac0 , 1.3*Rₕaux[1] )/(1 +_ten * _ten/coup.ten0) #velocity polimerization  
    _xₗ = xₗ*0.9 + control.Δt*vₗ #we introduce a slight relaxation for the membrane, decreases 5% x at the Boundary condition only
    xₗ = _xₗ 


    UXʷ,UVʷ,Vʷ,UXᵛ,UVᵛ,Xᵛ,Yᵛ,Xʳ,Yʳ,Uᵉ,Vᵉ,Vᴿ,Uᴿ,Ωᶜ,dΩᶜ,dΓ,nΓ,φ = update_all!(i,t,control.Δt,υₕ,msₕ)

    op_rho = AffineFEOperator(Arho,Brho,Uᴿ,Vᴿ)
    op_rac = AffineFEOperator(Arac,Brac,Uᴿ,Vᴿ)
    Rₕ = solve(op_rac)
    Rₕ_old = Rₕ
    ρₕ = solve(op_rho) 
    ρₕ_old = ρₕ 

    op_MCAb= AffineFEOperator(AMCAb,BMCAb,Uᴿ,Vᴿ)
    uh_MCAb=solve(op_MCAb)
    uh_MCAb_old=uh_MCAb
  
    op_MCAu= AffineFEOperator(AMCAu,BMCAu,Uᴿ,Vᴿ)
    uh_MCAu=solve(op_MCAu)
    uh_MCAu_old=uh_MCAu
     
    ract[i,:] = vcat(lazy_map(Rₕ,xΓ)...)[perm]
    rhot[i,:] = vcat(lazy_map(ρₕ,xΓ)...)[perm]
    MCAbt[i,:] = vcat(lazy_map(uh_MCAb,xΓ)...)[perm]
    σₐt[i,:] =    mech.sigmaₐ⁰ .+ mech.sigmaρ⁰ * rhot[i,:]  
    χt[i,:] =  mech.χ₀ .+ mech.χ*MCAbt[i,:] #χR(ract[i,:]) 
    _vt = υₕ⋅τ
    vt[i,:]  = vcat(lazy_map(_vt,xΓ)...)[perm] 
    _vnt = vcat(lazy_map(υₕ ⋅ nΓ,xΓ)...) 
    vnt[i,:] = _vnt[perm] 
    _xt = xₕ⋅τ
    xt[i,:]  = vcat(lazy_map(_xt,xΓ)...)[perm] 
    tent[i,:] = vcat(lazy_map(ten,xΓ)...)[perm]
    
    plotting("rac",ract[i,:],pPNG*"Rac_time/","$i")
    plotting("rho",rhot[i,:],pPNG*"Rho_time/","$i") 
    plotting("MCAb",MCAbt[i,:],pPNG*"MCAb_time/","$i") 
  end
  
  # plots_run_singlet(nΔt,vt*0.0006,xt*6,MCAbt,ract,rhot,pPNG,
  #  num_qpoints,π*R2*6,Δt₀*10000,control.T*10000,flat_alenΓ,σₐt,χt,tent*10)  #nΔt,vt,ract,rhot,pPNG,   partition,L,control.Δt,control.T,xplot
  plots_run_singlet(nΔt,vt,vnt,xt*6,MCAbt,ract,rhot,pPNG,
   num_qpoints,π*R2,control.Δt,control.T,flat_alenΓ,σₐt,χt,tent)  #nΔt,vt,ract,rhot,pPNG,   partition,L,control.Δt,control.T,xplot

end