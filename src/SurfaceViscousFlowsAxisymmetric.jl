include("Plots_RhoRacSinglet.jl")

function surface_viscous_flows_axisymmetric(
            dᵃ,α₀,Drac,αopto,dᵇ,β₀,Drho,βopto,wrac,
            domain::Tuple{Vararg{Float64}},
            ls::AlgoimCallLevelSetFunction,
            Pe::Float64,
            n::Int,
            Δt₀::Float64,
            T::Float64,
            σₐ⁰::Float64,
            χ₀::Float64,
            χ::Float64;
            initial_density::Function = verification,
            activity::Function = unit_activity_axisymmetric,
            order::Int = 2,
            γᶜ::Float64 = 1.0, 
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
    Yᵛ = MultiFieldFESpace([Vʷ,Vˡ,Vˡ])
    Xᵛ = MultiFieldFESpace([Uʷ,Uˡ,Uˡ])

    # Space to create homogeneous perturbation  
    # of constant concentration myosin field
    # TO-DO: To be deleted if needed
    Yʳ = MultiFieldFESpace([Vᵉ,Vˡ])
    Xʳ = MultiFieldFESpace([Uᵉ,Uˡ])

    Xᵛ,Yᵛ,Xʳ,Yʳ,Uᵉ,Vᵉ,Vᴿ,Uᴿ,dΩᶜ,dΓ,nΓ,φ

  end

 # Lets make output folders
  pVTU="./VTU/"*name
  mkpath(pVTU)
  pPNG="./PNG/"*name
  mkpath(pPNG)
  mkpath(pPNG*"Rac_time/") 
  mkpath(pPNG*"Rho_time/") 
  mkpath(pPNG*"Rac_time_initial/") 
  mkpath(pPNG*"Rho_time_initial/") 
 
  # Lets copy the code in the output folder to be able to check code used for each simulation
  cp(@__FILE__, pPNG*split(@__FILE__, "/")[end],force=true)

  # Time discretisation parameters
  t₀ = 0.0
  Δt = Δt₀
  u₀ = VectorValue(0.0,0.0)
  m₀ = 2.0
  R2=1
  R=1
  nΔt = trunc(Int,T/Δt+0.5)+1

  Xᵛ,Yᵛ,Xʳ,Yʳ,Uᵉ,Vᵉ,Vᴿ,Uᴿ,dΩᶜ,dΓ,nΓ,φ = update_all!(0,t₀,Δt,u₀,m₀)

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

  #Building the vectors used for introducing opto influence as an increase in α and β
  arclength(x) = R2 * atan(x[2],-x[1]) # Arc length for sphere
  α₀opto(x) = α₀
  β₀opto(x) = β₀

  α₀opto2(x) = α₀ + αopto * exp( -0.5 * ( arclength(x)-π*R2 )^2 / ((wrac)^2) )
  β₀opto2(x) = β₀ + βopto * exp( -0.5 * ( arclength(x))^2       / ((wrac)^2) )

  γ₀ = 0.1  / h # TODO: Eric reviews the scaling with h
  γ₀R =  0.1  / h 
  m₀opto(u,v) = ∫( u*v )dΓ
  s₀opto(u,v) = ∫( γ₀*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ
  s₀R(u,v) = ∫( γ₀R*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ

  A₀opto(u,v) = m₀opto(u,v) + s₀opto(u,v)
  bα₀opto(v) = m₀opto(α₀opto,v)
  bβ₀opto(v) = m₀opto(β₀opto,v)
  bα₀opto2(v) = m₀opto(α₀opto2,v)
  bβ₀opto2(v) = m₀opto(β₀opto2,v)

  op_α₀ = AffineFEOperator(A₀opto,bα₀opto,Uᴿ,Vᴿ)
  op_β₀ = AffineFEOperator(A₀opto,bβ₀opto,Uᴿ,Vᴿ)

  α₀v = solve(op_α₀)
  β₀v = solve(op_β₀)

  #Rac and Rho initialization
  uh_rac = interpolate_everywhere(0.0,Uᴿ) 
  uh_rho = interpolate_everywhere(4.0,Uᴿ) 
  uh_rac_old = uh_rac
  uh_rho_old = uh_rho
  a_rac, b_rac, a_rho, b_rho = rac_rho_weak_forms2(Δt,dᵃ,dᵇ,Drac,Drho,nΓ,dΓ)
 
  #SOLVE Rac AT t=0
  Arac(rac,w) = a_rac(rac,w) + s₀R(rac,w)
  Brac(w) = b_rac(w,uh_rho,α₀v,uh_rac_old) #(w,rho,α₀v,rac_old)
  op_rac = AffineFEOperator(Arac,Brac,Uᴿ,Vᴿ)
  uh_rac = solve(op_rac)
  uh_rac_old = uh_rac
  #SOLVE Rho AT t=0
  Arho(rho,w) = a_rho(rho,w) + s₀R(rho,w)
  Brho(w) = b_rho(w,uh_rac,β₀v,uh_rho_old) #(w,rac,β₀v,rho_old)
  op_rho = AffineFEOperator(Arho,Brho,Uᴿ,Vᴿ)
  uh_rho = solve(op_rho) 
  uh_rho_old = uh_rho
  
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
  vt = zeros(nΔt,num_qpoints)
  _vt = vcat(lazy_map(υₕ,xΓ)...) 
  _vt = _vt[perm]
  #vv = get_cell_dof_values(_vt)
  vt[1,:] .= √(_vt⋅_vt)

for ti in 1:100 
    op_rho = AffineFEOperator(Arho,Brho,Uᴿ,Vᴿ)
    uh_rho = solve(op_rho)
    uh_rho_old = uh_rho
    op_rac = AffineFEOperator(Arac,Brac,Uᴿ,Vᴿ)
    uh_rac = solve(op_rac)
    uh_rac_old = uh_rac
    sum_rac = ∑(∫(uh_rac)dΓ)
    sum_rho = ∑(∫(uh_rho)dΓ)
    ractt = vcat(lazy_map(uh_rac,xΓ)...)
    rhott = vcat(lazy_map(uh_rho,xΓ)...)
    ractt[:] = ractt[perm]
    rhott[:] = rhott[perm]
    plotting("rac",ractt[:],pPNG*"Rac_time_initial/","$ti")
    plotting("rho",rhott[:],pPNG*"Rho_time_initial/","$ti")
  end

  while t < T + tol
    if i > 100
      op_α₀ = AffineFEOperator(A₀opto,bα₀opto2,Uᴿ,Vᴿ)
      op_β₀ = AffineFEOperator(A₀opto,bβ₀opto2,Uᴿ,Vᴿ)
      α₀v = solve(op_α₀)
      β₀v = solve(op_β₀)
    end
    if i > 200
      op_α₀ = AffineFEOperator(A₀opto,bα₀opto,Uᴿ,Vᴿ)
      op_β₀ = AffineFEOperator(A₀opto,bβ₀opto,Uᴿ,Vᴿ)
      α₀v = solve(op_α₀)
      β₀v = solve(op_β₀)
    end
    
    @info "Time step $i, time $t and time step $Δt"

    assemᵛ = SparseMatrixAssembler(Tm,Tv,Xᵛ,Yᵛ)
    Aᵛ = nothing

    aᵛ,bᵛ = cortical_flow_problem_axisymmetric(
        uh_rho,uh_rac,dΩᶜ,dΓ,nΓ,γʷ,Pe,χ,χ₀,σₐ⁰)
    Aᵛ,Bᵛ = _assemble_problem(aᵛ,bᵛ,assemᵛ,Xᵛ,Yᵛ,Aᵛ)
    υₕ,_ = _solve_problem(Aᵛ,Bᵛ,Xᵛ,ps)

    writesol && postprocess_all(φ,dΩᶜ.quad.trian,
      uh_rac,uh_rho,υₕ,i=i,of=output_frequency,name=pVTU)

    msₕ = get_maximum_magnitude_with_dirichlet(υₕ)

    i = i + 1
    t = t + Δt

    Xᵛ,Yᵛ,Xʳ,Yʳ,Uᵉ,Vᵉ,Vᴿ,Uᴿ,dΩᶜ,dΓ,nΓ,φ = update_all!(i,t,Δt,υₕ,msₕ)

    op_rho = AffineFEOperator(Arho,Brho,Uᴿ,Vᴿ)
    op_rac = AffineFEOperator(Arac,Brac,Uᴿ,Vᴿ)
    uh_rac = solve(op_rac)
    uh_rac_old = uh_rac
    uh_rho = solve(op_rho) 
    uh_rho_old = uh_rho

    ract[i,:] = vcat(lazy_map(uh_rac,xΓ)...)
    rhot[i,:] = vcat(lazy_map(uh_rho,xΓ)...) 
    ract[i,:] = ract[i,perm]
    rhot[i,:] = rhot[i,perm]

    # fa(μ, ρ,R) = ∫( ( -(divᶜ(μ,nΓ)+μ⋅iy)*(∇ᵈ(ρ,nΓ)⋅(TensorValue(0.0,-1.0,1.0,0.0)⋅nΓ) - ∇ᵈ(R,nΓ)⋅(TensorValue(0.0,-1.0,1.0,0.0)⋅nΓ)) ) * χ₀ )dΓ
    # papa(v) = fa(v,uh_rho,uh_rac)
    # aux_ = AffineFEOperator(A₀opto,papa,Uᴿ,Vᴿ) 
    # aux = solve(aux_)
    # plotting("gradients",get_cell_dof_values(aux),pPNG,"$i")
    plotting("rac",ract[i,:],pPNG*"Rac_time/","$i")
    plotting("rho",rhot[i,:],pPNG*"Rho_time/","$i")
  end
  
  plots_run_singlet(nΔt,vt,ract,rhot,pPNG,
   num_qpoints,π*R2,Δt₀,T,flat_alenΓ)  #nΔt,vt,ract,rhot,pPNG,   partition,L,Δt,T,xplot

end




function surface_viscous_flows_axisymmetric_conserved(
            dᵃ,α₀,Drac,αopto,dᵇ,β₀,Drho,βopto,wrac,
            domain::Tuple{Vararg{Float64}},
            ls::AlgoimCallLevelSetFunction,
            Pe::Float64,
            n::Int,
            Δt₀::Float64,
            T::Float64,
            rac_total::Float64, rho_total::Float64,
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
    Yᵛ = MultiFieldFESpace([Vʷ,Vˡ,Vˡ])
    Xᵛ = MultiFieldFESpace([Uʷ,Uˡ,Uˡ])

    # Space to create homogeneous perturbation  
    # of constant concentration myosin field
    # TO-DO: To be deleted if needed
    Yʳ = MultiFieldFESpace([Vᵉ,Vˡ])
    Xʳ = MultiFieldFESpace([Uᵉ,Uˡ])

    Xᵛ,Yᵛ,Xʳ,Yʳ,Uᵉ,Vᵉ,Vᴿ,Uᴿ,dΩᶜ,dΓ,nΓ,φ

  end

 # Lets make output folders
  pVTU="./VTU/"*name
  mkpath(pVTU)
  pPNG="./PNG/"*name
  mkpath(pPNG)
  mkpath(pPNG*"Rac_time/") 
  mkpath(pPNG*"Rho_time/") 
  mkpath(pPNG*"Rac_time_initial/") 
  mkpath(pPNG*"Rho_time_initial/") 
 
  # Lets copy the code in the output folder to be able to check code used for each simulation
  cp(@__FILE__, pPNG*split(@__FILE__, "/")[end],force=true)

  # Time discretisation parameters
  t₀ = 0.0
  Δt = Δt₀
  u₀ = VectorValue(0.0,0.0)
  m₀ = 2.0
  R2=1
  R=1
  nΔt = trunc(Int,T/Δt+0.5)+1

  Xᵛ,Yᵛ,Xʳ,Yʳ,Uᵉ,Vᵉ,Vᴿ,Uᴿ,dΩᶜ,dΓ,nΓ,φ = update_all!(0,t₀,Δt,u₀,m₀)

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

  #Building the vectors used for introducing opto influence as an increase in α and β
  arclength(x) = R2 * atan(x[2],-x[1]) # Arc length for sphere
  α₀opto(x) = α₀
  β₀opto(x) = β₀

  α₀opto2(x) = α₀ + αopto * exp( -0.5 * ( arclength(x)-π*R2 )^2 / ((wrac)^2) )
  β₀opto2(x) = β₀ + βopto * exp( -0.5 * ( arclength(x))^2       / ((wrac)^2) )

  γ₀ = 0.1  / h # TODO: Eric reviews the scaling with h
  γ₀R =  0.1  / h 
  m₀opto(u,v) = ∫( u*v )dΓ
  s₀opto(u,v) = ∫( γ₀*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ
  s₀R(u,v) = ∫( γ₀R*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ

  A₀opto(u,v) = m₀opto(u,v) + s₀opto(u,v)
  bα₀opto(v) = m₀opto(α₀opto,v)
  bβ₀opto(v) = m₀opto(β₀opto,v)
  bα₀opto2(v) = m₀opto(α₀opto2,v)
  bβ₀opto2(v) = m₀opto(β₀opto2,v)

  op_α₀ = AffineFEOperator(A₀opto,bα₀opto,Uᴿ,Vᴿ)
  op_β₀ = AffineFEOperator(A₀opto,bβ₀opto,Uᴿ,Vᴿ)

  α₀v = solve(op_α₀)
  β₀v = solve(op_β₀)

  #Rac and Rho initialization
  uh_rac = interpolate_everywhere(0.0,Uᴿ) 
  uh_rho = interpolate_everywhere(0.0,Uᴿ) 
  uh_rac_i = interpolate_everywhere(rac_total/(π*R2),Uᴿ) 
  uh_rho_i = interpolate_everywhere(rho_total/(π*R2),Uᴿ) 
  uh_rac_old = uh_rac
  uh_rho_old = uh_rho
  a_rac, b_rac, a_rho, b_rho, a_rac_i, b_rac_i, a_rho_i, b_rho_i = rac_rho_weak_forms_conserved(Δt,dᵃ,dᵇ,Drac,Drho,nΓ,dΓ)
 
  #SOLVE Rac AT t=0
  Arac(rac,w) = a_rac(rac,w) + s₀R(rac,w)
  Brac(w) = b_rac(w,uh_rho,α₀v,uh_rac_i,uh_rac_old) #(w,rho,α₀v,rac_old)
  op_rac = AffineFEOperator(Arac,Brac,Uᴿ,Vᴿ)
  uh_rac = solve(op_rac)
  uh_rac_old = uh_rac
  sum_rac = ∑(∫(uh_rac)dΓ)
  #Solve Rac inactive at t=0
  Arac_i(rac_i,w) = a_rac_i(rac_i,w) + s₀R(rac_i,w)
  Brac_i(w) = b_rac_i(w,sum_rac) #(w,rho,α₀v,rac_old)
  op_rac_i = AffineFEOperator(Arac_i,Brac_i,Uᴿ,Vᴿ)
  uh_rac_i = solve(op_rac_i)
  #SOLVE Rho AT t=0
  Arho(rho,w) = a_rho(rho,w) + s₀R(rho,w)
  Brho(w) = b_rho(w,uh_rac,β₀v,uh_rho_old) #(w,rac,β₀v,rho_old)
  op_rho = AffineFEOperator(Arho,Brho,Uᴿ,Vᴿ)
  uh_rho = solve(op_rho) 
  uh_rho_old = uh_rho
  sum_rho = ∑(∫(uh_rho)dΓ)
  #Solve Rho inactive at t=0
  Arho_i(rho_i,w) = a_rho_i(rho_i,w) + s₀R(rho_i,w)
  Brho_i(w) = b_rho_i(w,sum_rho) #(w,rho,α₀v,rac_old)
  op_rho_i = AffineFEOperator(Arho_i,Brho_i,Uᴿ,Vᴿ)
  uh_rho_i = solve(op_rho_i)
  
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
  rac_it = zeros(nΔt,num_qpoints)
  rho_it = zeros(nΔt,num_qpoints)
  vt = zeros(nΔt,num_qpoints)
  _vt = vcat(lazy_map(υₕ,xΓ)...) 
  _vt = _vt[perm]
  #vv = get_cell_dof_values(_vt)
  vt[1,:] .= √(_vt⋅_vt)

for ti in 1:100 
    op_rho = AffineFEOperator(Arho,Brho,Uᴿ,Vᴿ)
    uh_rho = solve(op_rho)
    uh_rho_old = uh_rho
    op_rac = AffineFEOperator(Arac,Brac,Uᴿ,Vᴿ)
    uh_rac = solve(op_rac)
    uh_rac_old = uh_rac

    op_rho_i = AffineFEOperator(Arho_i,Brho_i,Uᴿ,Vᴿ)
    uh_rho_i = solve(op_rho_i) 
    op_rac_i = AffineFEOperator(Arac_i,Brac_i,Uᴿ,Vᴿ)
    uh_rac_i = solve(op_rac_i) 

    sum_rac = ∑(∫(uh_rac)dΓ)
    sum_rho = ∑(∫(uh_rho)dΓ)
    ractt = vcat(lazy_map(uh_rac,xΓ)...)
    rhott = vcat(lazy_map(uh_rho,xΓ)...)
    ractt[:] = ractt[perm]
    rhott[:] = rhott[perm]
    plotting("rac",ractt[:],pPNG*"Rac_time_initial/","$ti")
    plotting("rho",rhott[:],pPNG*"Rho_time_initial/","$ti")
  end

  while t < T + tol
    if i > 100
      op_α₀ = AffineFEOperator(A₀opto,bα₀opto2,Uᴿ,Vᴿ)
      op_β₀ = AffineFEOperator(A₀opto,bβ₀opto2,Uᴿ,Vᴿ)
      α₀v = solve(op_α₀)
      β₀v = solve(op_β₀)
    end
    if i > 200
      op_α₀ = AffineFEOperator(A₀opto,bα₀opto,Uᴿ,Vᴿ)
      op_β₀ = AffineFEOperator(A₀opto,bβ₀opto,Uᴿ,Vᴿ)
      α₀v = solve(op_α₀)
      β₀v = solve(op_β₀)
    end
    
    @info "Time step $i, time $t and time step $Δt"

    assemᵛ = SparseMatrixAssembler(Tm,Tv,Xᵛ,Yᵛ)
    Aᵛ = nothing

    aᵛ,bᵛ = cortical_flow_problem_axisymmetric(
        uh_rho,dΩᶜ,dΓ,nΓ,γʷ,Pe,χ,σₐ⁰)
    Aᵛ,Bᵛ = _assemble_problem(aᵛ,bᵛ,assemᵛ,Xᵛ,Yᵛ,Aᵛ)
    υₕ,_ = _solve_problem(Aᵛ,Bᵛ,Xᵛ,ps)

    writesol && postprocess_all(φ,dΩᶜ.quad.trian,
      eₕ,υₕ,i=i,of=output_frequency,name=name)

    msₕ = get_maximum_magnitude_with_dirichlet(υₕ)

    i = i + 1
    t = t + Δt

    Xᵛ,Yᵛ,Xʳ,Yʳ,Uᵉ,Vᵉ,Vᴿ,Uᴿ,dΩᶜ,dΓ,nΓ,φ = update_all!(i,t,Δt,υₕ,msₕ)

    assemᵉ = SparseMatrixAssembler(Tm,Tv,Uᵉ,Vᵉ)
    aᵉ,bᵉ = transport_problem_axisymmetric(
      υₕ,eₕ,dΓ,dΩᶜ,nΓ,Δt,γᵉ,τᵈkₒ)
    opᵉ = AffineFEOperator(aᵉ,bᵉ,Uᵉ,Vᵉ,assemᵉ)
    eₕ = solve(ps,opᵉ)

    op_rho = AffineFEOperator(Arho,Brho,Uᴿ,Vᴿ)
    op_rac = AffineFEOperator(Arac,Brac,Uᴿ,Vᴿ)
    uh_rac = solve(op_rac)
    uh_rac_old = uh_rac
    uh_rho = solve(op_rho) 
    uh_rho_old = uh_rho
    #vt[i+1,:] = vcat(lazy_map(υₕ,xΓ)...) 
    #vt[i+1,:] = vt[i+1,perm]

    op_rho_i = AffineFEOperator(Arho_i,Brho_i,Uᴿ,Vᴿ)
    uh_rho_i = solve(op_rho_i) 
    op_rac_i = AffineFEOperator(Arac_i,Brac_i,Uᴿ,Vᴿ)
    uh_rac_i = solve(op_rac_i) 

    ract[i,:] = vcat(lazy_map(uh_rac,xΓ)...)
    rhot[i,:] = vcat(lazy_map(uh_rho,xΓ)...) 
    ract[i,:] = ract[i,perm]
    rhot[i,:] = rhot[i,perm]

    plotting("rac",ract[i,:],pPNG*"Rac_time/","$i")
    plotting("rho",rhot[i,:],pPNG*"Rho_time/","$i")
  end
  
  plots_run_singlet(nΔt,vt,ract,rhot,pPNG,
   num_qpoints,π*R2,Δt₀,T,flat_alenΓ)  #nΔt,vt,ract,rhot,pPNG,   partition,L,Δt,T,xplot

end