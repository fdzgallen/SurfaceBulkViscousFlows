# Exemple on fixed surface with variable thickness

# TO-DO: Pass from 2D to 2D axisymmetric
function variable_thickness_cortical_flows(tf::TraceFEM,
                                           phi::LevelSetManifold,
                                           prob::VariableThickness,
                                           n::Int,
                                           order::Int)

  domain = phi.domain
  partition = (n,n)
  h = (domain[2]-domain[1])/n
  bgmodel = CartesianDiscreteModel(domain,partition)
  Ωbg = Triangulation(bgmodel)

  # Algoim quadrature on the background model
  degree = order == 1 ? 3 : 2*order
  levelset = phi.ls
  squad = Quadrature(algoim,levelset,degree)
  dΓbg = Measure(Ωbg,squad,data_domain_style=PhysicalDomain())

  # Active triangulation and its associated standard quadrature
  Ωc = Triangulation(Ωbg,is_cell_active(dΓbg))
  dΩc = Measure(Ωc,2*order)

  # Algoim quadrature and normal on embedded manifold
  dΓ = restrict_measure(dΓbg,Ωc)
  n_Γ = normal(levelset,Ωc)

  # Reference FEs
  D = num_dims(bgmodel)
  reffeᵘ = ReferenceFE(lagrangian,VectorValue{D,Float64},order)
  reffeᵖ = ReferenceFE(lagrangian,Float64,order-1)
  reffeᵉ = ReferenceFE(lagrangian,Float64,order-1)

  # Test FE spaces
  V = TestFESpace(Ωc,reffeᵘ)
  _Q = TestFESpace(Ωc,reffeᵖ)
  Q = ZeroMeanFESpace(_Q,dΓ)
  F = TestFESpace(Ωc,reffeᵉ) # For thickness

  # Trial FE spaces
  U = TrialFESpace(V)
  P = TrialFESpace(Q)
  E = TransientTrialFESpace(F) # For thickness (no need to be transient)

  # Multifield FE spaces
  Y = MultiFieldFESpace([V,Q])
  X = TransientMultiFieldFESpace([U,P])
  # Y = MultiFieldFESpace([V,Q,F])
  # X = TransientMultiFieldFESpace([U,P,E])

  u,p,e,fᵘ,fᵖ,fᵉ = 
    cortical_flow_problem_functions(phi,phi.name,prob)

  # TO-DO: Inspect weak forms because they hint how to
  # adapt the surface_viscous_flows_axisymmetric to variable thickness

  # *** CORTICAL FLOW ***
  # ** grad-grad term **
  # RMK. Adding zeroth term to circumvent nontrivial kernel
  μ = phi.μ; σᵛ(ε) = 2*μ*ε; σᵈ(v,n) = σᵛ∘εᵈ(v,n)
  # Key Takeaway 1 -> Everything premultiplied by thickness e in cortical flows
  a(u,v,e) = ∫( e * ( σᵈ(u,n_Γ)⊙∇ᵈ(v,n_Γ) 
                 - 2*(σᵈ(u,n_Γ)⋅n_Γ)⋅(εᵈ(v,n_Γ)⋅n_Γ) + u⋅v ) )dΓ
  # -> Do not pay attention to what follows
  # ** pressure term **
  b(v,q,e) = ∫( ( v⋅(e*∇ᵈ(q,n_Γ)+q*∇ᵈ(e,n_Γ)) ) )dΓ
  # b(v,q,e) = ∫( ( tr(∇ᵈ(v,n_Γ))*e*q ) )dΓ
  # ** tangent penalty **
  η = tf.stab_coeff/(h^2.0)
  k(u,v,e) = ∫( η*e*((u⋅n_Γ)*(v⋅n_Γ)) )dΓ
  # ** u-stabilisation **
  γᵘ = tf.stab_coeff/h
  sᵘ(u,v,e) = ∫( γᵘ*e*((n_Γ⋅∇(u))⊙(n_Γ⋅∇(v))) )dΩc
  # ** p-stabilisation **
  γᵖ = tf.stab_coeff*h
  sᵖ(p,q,e) = ∫( γᵖ*e*((n_Γ⋅∇(p))⊙(n_Γ⋅∇(q))) )dΩc

  # *** THICKNESS EVOLUTION ***
  # TO-DO: Do not implement with time derivative
  # implement analogous to transport problem in
  # surface_viscous_flows_axisymmetric.jl
  # ** time derivative term **
  mₜ(e,ε) = ∫( ∂t(e)*ε )dΓ
  # ** mass term **
  m(t,e,ε) = ∫( e*ε )dΓ
  # ** convective term **
  c(e,ε,u) =  ∫( 0.5 * ( 
    (u⋅∇ᵈ(e,n_Γ))*ε - e*(u⋅∇ᵈ(ε,n_Γ)) + tr(∇ᵈ(u,n_Γ))*e*ε ) )dΓ
  # ** reaction term **
  kᵈ = prob.kᵈ
  r(e,ε) = ∫( kᵈ*e*ε )dΓ
  # ** e-stabilisation **
  γᵉ = tf.stab_coeff*h # [?] Scale with velocities?
  sᵉ(e,ε) = ! tf.ghost_penalty ? ∫( γᵉ*((n_Γ⋅∇(e))⊙(n_Γ⋅∇(ε))) )dΩc : begin
    Λ = SkeletonTriangulation(Ωc)
    dΛ = Measure(Λ,2*order)
    n_Λ = get_normal_vector(Λ)
    ∫( γᵉ*jump(n_Λ⋅∇(e))*jump(n_Λ⋅∇(ε)) )dΛ 
  end
  # ** source term **
  vᵖ = prob.vᵖ
  l(t,v,q,e,ε) = ∫( vᵖ*ε )dΓ + 
    ∫( e*(fᵘ(t)⊙v) )dΓ + ∫( e*(fᵖ(t)*q) )dΓ + ∫( fᵉ(t)*ε )dΓ

  _aᵘ(u,v,p,q,e) = a(u,v,e) + sᵘ(u,v,e) + k(u,v,e) + 
                   b(v,p,e) - b(u,q,e) + sᵖ(p,q,e)
  _bᵘ(t,v,q) = ∫( fᵘ(t)⊙v )dΓ + ∫( fᵖ(t)*q )dΓ

  # This is to integrate thickness evolution in time
  # using a standard ODE solver, but we don't do it
  # like this in the surface_viscous_flows_axisymmetric.jl
  uₕ = interpolate_everywhere(u,U)
  aᵉ(t,e,ε) = r(e,ε) + c(e,ε,uₕ) + sᵉ(e,ε) # Send flow to thickness problem
  bᵉ(t,ε) = ∫( vᵖ*ε )dΓ + ∫( fᵉ(t)*ε )dΓ
  opᵉ = TransientAffineFEOperator(m,aᵉ,bᵉ,E,F)

  ls = LUSolver()
  Δt = 0.1
  θ = 1.0
  ode_solver = ThetaMethod(ls,Δt,θ)

  t₀ = 0.0; T = 1.0
  e₀ = interpolate_everywhere(e(t₀),E(t₀))
  eₜ = solve(ode_solver,opᵉ,e₀,t₀,T)

  l2(w)   = ∑( ∫( w⊙w )dΓ )
  h1(w)   = ∑( ∫( ∇ᵈ(w,n_Γ)⊙∇ᵈ(w,n_Γ) )dΓ )
  mass(w) = ∑( ∫( w )dΓ )

  el2ₑ = 0.0; eh1ₑ = 0.0;
  el2ᵤ = 0.0; eh1ᵤ = 0.0;
  el2ₚ = 0.0; eh1ₚ = 0.0;

  # x = dΓ.quad.cell_point.values
  # x = lazy_map(Reindex(x),dΓ.quad.cell_point.ptrs)

  # Key takeaway 2 -> Adapt this idea to your code
  # To solve the viscous flow problem with variable thickness
  # Solve in a decoupled way the thickness and the flow
  # And send to each one the FE functions of the other variables
  # Staggered scheme
  # (1) Solve for thickness e
  # (2) Solve for velocity u and pressure p
  for (i,(eₕ,t)) in enumerate(eₜ)
    eₑ   = e(t) - eₕ
    el2ₑ = el2ₑ + l2(eₑ)
    eh1ₑ = eh1ₑ + h1(eₑ)
    aᵘ((u,p),(v,q)) = _aᵘ(u,v,p,q,eₕ) # Send thickness to flow problem
    bᵘ((v,q)) = _bᵘ(t,v,q)
    opᵘ  = AffineFEOperator(aᵘ,bᵘ,X,Y)
    uₕ, pₕ = solve(opᵘ)
    eᵤ   = u - uₕ
    el2ᵤ = el2ᵤ + l2(eᵤ)
    eh1ᵤ = eh1ᵤ + h1(eᵤ)
    eₚ   = p - pₕ
    el2ₚ = el2ₚ + l2(eₚ)
    eh1ₚ = eh1ₚ + h1(eₚ)
    # ux = lazy_map(uₕ,x)
    # px = lazy_map(pₕ,x)
    # eux = lazy_map(eᵤ,x)
    # epx = lazy_map(eₚ,x)
    # writevtk(x,"sres_$i",nodaldata=["uₕ"=>ux,"euₕ"=>eux,"pₕ"=>px,"epₕ"=>epx])
  end

  el2ₑ = √(Δt*el2ₑ); eh1ₑ = √(Δt*eh1ₑ)
  el2ᵤ = √(Δt*el2ᵤ); eh1ᵤ = √(Δt*eh1ᵤ)
  el2ₚ = √(Δt*el2ₚ); eh1ₚ = √(Δt*eh1ₚ)
  @show el2ᵤ, eh1ᵤ, el2ₚ, eh1ₚ, el2ₑ, eh1ₑ

end
