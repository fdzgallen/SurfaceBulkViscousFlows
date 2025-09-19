function surface_viscous_flows_axisymmetric(
            domain::Tuple{Vararg{Float64}},
            ls::AlgoimCallLevelSetFunction,
            Pe::Float64,
            n::Int,
            Δt₀::Float64,
            T::Float64;
            initial_density::Function = verification,
            activity::Function = unit_activity_axisymmetric,
            order::Int = 2,
            γᶜ::Float64 = 1.0,
            χ::Float64 = 1.0,
            τᵈkₒ::Float64 = 10.0,
            writesol::Bool = true,
            output_frequency::Int = 1,
            redistance_frequency::Int = 1,
            name::String = "plt")

  # Background geometry
  cells = (n,div(n,2))
  h = (domain[2]-domain[1])/n
  bgmodel = CartesianDiscreteModel(domain,cells)
  Ω = Triangulation(bgmodel)

  # Buffer of active model and integration objects
  degree = order < 3 ? 3 : 2*order
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
        Vbg = TestFESpace(Ω,ReferenceFE(lagrangian,Float64,order))
        _φ₋ = interpolate_everywhere(ls.φ,Vbg)
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
        _φ₋ = compute_distance_fe_function(bgmodel,Vbg,φ₋,order,cppdegree=3)
        φ₋  = AlgoimCallLevelSetFunction(_φ₋,∇(_φ₋))
      end

      cp₋ = compute_closest_point_projections(Vbg,φ₋,order,
              cppdegree=3,trim=true,limitstol=1.0e-2)

      # Current time surface and bulk measures
      squad = Quadrature(algoim,φ₋,degree,phase=CUT)
      s_cell_quad,is_c₋ = CellQuadratureAndActiveMask(bgmodel,squad)

      # Surface narrow-band triangulation
      δ₋ = 2.0 * mv₋₂ * dt
      _,is_nᶜ = narrow_band_triangulation(Ω,_φ₋,Vbg,is_c₋,δ₋)

      Ωᶜ,dΓ = TriangulationAndMeasure(Ω,s_cell_quad,is_nᶜ,is_c₋)

      dΩᶜ = Measure(Ωᶜ,2*order)
      nΓ = normal(φ₋,Ω)

      # Update buffer
      buffer[] = ( Ωᶜ=Ωᶜ,dΩᶜ=dΩᶜ,dΓ=dΓ,nΓ=nΓ,cp₋=cp₋,φ₋=φ₋,t=t,Vbg=Vbg )
      return true

    end

  end

  # Reference FEs
  N = num_dims(bgmodel)
  reffeʷ = ReferenceFE(lagrangian,VectorValue{N,Float64},order-1)
  reffeᵉ = ReferenceFE(lagrangian,Float64,order-1)

  function update_all!(i::Int,t::Real,dt::Real,disp,val::Real)

    update_buffer!(i,t,dt,disp,val)

    # Triangulations and aggregates
    Ωᶜ = buffer[].Ωᶜ

    # Measures and normal
    dΩᶜ = buffer[].dΩᶜ
    dΓ  = buffer[].dΓ
    nΓ  = buffer[].nΓ
    φ   = buffer[].φ₋

    # Test FE spaces

    # u-surface
    Vʷ = TestFESpace(Ωᶜ,reffeʷ,dirichlet_tags=[5],
                               dirichlet_masks=[(false,true)])
    # e-surface
    Vᵉ = TestFESpace(Ωᶜ,reffeᵉ)
    # Lagrange multipliers
    Vˡ = ConstantFESpace(bgmodel)

    # Trial FE spaces
    Uʷ = TrialFESpace(Vʷ)
    Uᵉ = TrialFESpace(Vᵉ)
    Uˡ = TrialFESpace(Vˡ)

    # Multifield FE spaces
    Yᵛ = MultiFieldFESpace([Vʷ,Vˡ,Vˡ])
    Xᵛ = MultiFieldFESpace([Uʷ,Uˡ,Uˡ])

    # Space to create homogeneous perturbation  
    # of constant concentration myosin field
    # TO-DO: To be deleted if needed
    Yʳ = MultiFieldFESpace([Vᵉ,Vˡ])
    Xʳ = MultiFieldFESpace([Uᵉ,Uˡ])

    Xᵛ,Yᵛ,Xʳ,Yʳ,Uᵉ,Vᵉ,dΩᶜ,dΓ,nΓ,φ

  end

  # Time discretisation parameters
  t₀ = 0.0
  Δt = Δt₀
  u₀ = VectorValue(0.0,0.0)
  m₀ = 2.0

  Xᵛ,Yᵛ,Xʳ,Yʳ,Uᵉ,Vᵉ,dΩᶜ,dΓ,nΓ,φ = update_all!(0,t₀,Δt,u₀,m₀)

  # *** WEAK FORM PARAMETERS ***
  ξ(e) = 2.0 * e*e / ( 1.0 + e*e )
  # ** u-stabilisation **
  γʷ = γᶜ/h
  # ** e-stabilisation **
  γᵉ = γᶜ/h

  # Compute initial condition for myosin
  _eₕ = initial_density(Uᵉ,Xʳ,Yʳ,dΓ,dΩᶜ,nΓ)
  eₕ = interpolate_everywhere(_eₕ,Uᵉ)

  # Compute initial condition for surface velocity
  _υₕ(x) = VectorValue(0.0,0.0)
  υₕ  = interpolate_everywhere(_υₕ,Xᵛ[1])

  Tm = SparseMatrixCSR{0,PetscScalar,PetscInt}
  Tv = Vector{PetscScalar}
  ps = PETScLinearSolver(mykspsetup)

  i = 0
  t = t₀
  
  tol = 1e-8

  while t < T + tol

    @info "Time step $i, time $t and time step $Δt"

    assemᵛ = SparseMatrixAssembler(Tm,Tv,Xᵛ,Yᵛ)
    Aᵛ = nothing

    aᵛ,bᵛ = cortical_flow_problem_axisymmetric(
        eₕ,dΩᶜ,dΓ,nΓ,γʷ,Pe,χ,activity)
    Aᵛ,Bᵛ = _assemble_problem(aᵛ,bᵛ,assemᵛ,Xᵛ,Yᵛ,Aᵛ)
    υₕ,_ = _solve_problem(Aᵛ,Bᵛ,Xᵛ,ps)

    writesol && postprocess_all(φ,dΩᶜ.quad.trian,
      eₕ,υₕ,i=i,of=output_frequency,name=name)

    msₕ = get_maximum_magnitude_with_dirichlet(υₕ)

    i = i + 1
    t = t + Δt

    Xᵛ,Yᵛ,Xʳ,Yʳ,Uᵉ,Vᵉ,dΩᶜ,dΓ,nΓ,φ = update_all!(i,t,Δt,υₕ,msₕ)

    assemᵉ = SparseMatrixAssembler(Tm,Tv,Uᵉ,Vᵉ)
    aᵉ,bᵉ = transport_problem_axisymmetric(
      υₕ,eₕ,dΓ,dΩᶜ,nΓ,Δt,γᵉ,τᵈkₒ)
    opᵉ = AffineFEOperator(aᵉ,bᵉ,Uᵉ,Vᵉ,assemᵉ)
    eₕ = solve(ps,opᵉ)

  end

end