# Standard code to run a simulation using the mechanochemical model presented in article "INSERT DOI"
# The system will start with flat Rho and Rac at an steady state
# Different mechanical parametes can eb changed in the main section of the code.
# One can work with pure local inhibition by changing the following mechanical parameters:
# vCTE=0   α=0   β=0   σₐ₀=0
# Written by Andreu F Gallen working in Turlier lab and in collaboration with Orion Weiner's lab

include("Plots_RhoRacA.jl")

function plotting(ylab,po,pPNG,i)
  plot(po)
  xlabel!("ξ[μm]")
  ylabel!(ylab)
  savefig(pPNG*"rho"*i*".png")
end

function conservation(sezrin_u,sezrin_b,Minitial)
  return 2.0*(Minitial-(sezrin_u+sezrin_b))
end
 
function threshold(x,x₀,xth)
  return  (0.5 * (tanh∘(x/x₀ - xth/x₀)+1))
end
function threshold2(x,x₀,xth)
  return  (0.5 * (tanh.(x/x₀ .- xth/x₀).+1))
end

function MCA_bound_unbound_weak_forms(Δt,kon,koff,R2,D,nΓ,dΓ)
    # Now for MAC bound
  # mass term for the temporal evolution MCA_b
  mMCA(Δt,MCA_b,w) = ∫( ( (MCA_b*w)/Δt )*y )dΓ
  aezrin_b(MCA_b,v,w) = ∫( ( (MCA_b*w)/Δt )*y )dΓ    + ∫( ( 0.01 * (∇ᵈ(MCA_b,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ + 
    ∫( ( koff * ( MCA_b * w ) )*y )dΓ +
    ∫( ( w * ( v * (∇ᵈ(MCA_b,nΓ)⋅(TensorValue(0.0,-1.0,1.0,0.0)⋅nΓ)) + # ∫( ( w * ( v * (∇ᵈ(MCA_b,nΓ)⋅(VectorValue(1.0,1.0))) + #
    MCA_b * (∇ᵈ(v,nΓ)⋅(TensorValue(0.0,-1.0,1.0,0.0)⋅nΓ)) ) )*y )dΓ #+ #MCA_b * (∇ᵈ(v,nΓ)⋅(VectorValue(1.0,1.0))) ) )*y )dΓ + #
    #∫( ( λᵇ * ( ( MCA_b*MCA_b*MCA_b ) * w ) )*y )dΓ
  bezrin_b(w,MCA_u,ezrin_b_old,λ) = ∫( ( kon * ( MCA_u * w ) + (ezrin_b_old*w)/Δt )*y )dΓ  +
    ∫( ( λ/(π*R2) * ( kon / (kon+koff) )* w )*y )dΓ # mMCA(Δt,ezrin_b_old,w)  

  #Now for MCA unbound
  # mass term for the temporal evolution MCA_b0
  aezrin_u(MCA_u,w) = ∫( ( (MCA_u*w)/Δt )*y )dΓ + ∫( ( D * ( ∇ᵈ(MCA_u,nΓ)⋅∇ᵈ(w,nΓ) ) )*y )dΓ + 
    ∫( ( kon * ( MCA_u * w ) )*y )dΓ #+ 
   # MCA_u * ( (∇ᵈ(x,nΓ)⋅(TensorValue(0.0,-1.0,1.0,0.0)⋅nΓ)) - (∇ᵈ(x_old,nΓ)⋅(TensorValue(0.0,-1.0,1.0,0.0)⋅nΓ)) ) / Δt ) )*y )dΓ # + 
    #∫( ( λᵇ * ( ( MCA_u*MCA_u*MCA_u ) * w ) )*y )dΓ
  bezrin_u(w,MCA_b,ezrin_u_old,λ) = ∫( ( koff * (MCA_b*w))*y )dΓ + ∫( ( (ezrin_u_old*w)/Δt )*y )dΓ  + 
    ∫( ( λ/(π*R2) * (koff/(kon+koff))* w )*y )dΓ # mMCA(Δt,ezrin_u_old,w) 

  mMCA,aezrin_b,bezrin_b,aezrin_u,bezrin_u
end

function rac_rho_weak_forms(Δt,dᵃ,dᵇ,Drac,Drho,α,β,α₀v,β₀v,rho0,ezrin_bth,sig0,tenth,a_t,b_t,uh_rac_old,uh_rho_old,nΓ,dΓ,mMCA)
  #DEFINING the equations for Rac and Rho
  a_rac(rac,w,rho,MCA_b,α₀v) = (1/dᵃ) * mMCA(Δt,rac,w) + 
    ∫( ( Drac * (∇ᵈ(rac,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ + ∫( ( w*rac )*y )dΓ #+ ∫( ( w*rac * ( α₀v/(1+rho*rho) + α*((0.5 - threshold(MCA_b,rho0,ezrin_bth))) / (1+rho*rho) ) )*y )dΓ
  b_rac(w,rho,MCA_b,α₀v,rac_old) = (1/dᵃ) * mMCA(Δt,rac_old,w) + 
    ∫( ( w*(α₀v)/(1+rho*rho) )*y )dΓ  #∫( ( w*(a_t)*(α₀v/(1+rho*rho) + α*((0.5 - threshold(MCA_b,rho0,ezrin_bth)))/(1+rho*rho)) )*y )dΓ  
 
  a_rho(rho,w,rac,β₀v) = (1/dᵇ)*mMCA(Δt,rho,w) + 
    ∫( ( Drho * (∇ᵈ(rho,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ + 
    ∫( ( w*rho )*y )dΓ# +  ∫( ( w*rho * ( β₀v/(1+rac*rac)))*y )dΓ + ∫( (w*rho*(β*threshold(ten,sig0,tenth)/(1+rac*rac) ) )*y )dΓ # +
    # ∫( ( λʳᴬ*((rho*rho*rho)*w) )*y )dΓ  
  b_rho(w,rac,β₀v,rho_old) = (1/dᵇ)*mMCA(Δt,rho_old,w) + ∫( ( w*(β₀v)/(1+rac*rac) )*y )dΓ #+ ∫( ( w*b_t*(β₀v/(1+rac*rac)) + β*threshold(ten,sig0,tenth)/(1+rac*rac) )*y )dΓ

  a_rac, b_rac, a_rho, b_rho
end

#function to run a single simulation with a few given parameters
function run_singlet_axisymmetric(χ,η,T,Δt,part,
    L,simulation,wrac,αopto,βopto,kon,koff,M0,α₀,β₀,
    σₐ₀,D,Drac,Drho,rho0,ten0,a_t,b_t,α,β,dᵃ,dᵇ,
    sig0,tenth,ezrin_bth,topto,R,R2,L2)

  # Time discretisation parametersi
  t₀  = 0.0
  t   = t₀  
  nΔt = trunc(Int,T/Δt)

  # Level set function implicitly describing a unit sphere
  #R = 10.0 # Defined in main now
  φ = AlgoimCallLevelSetFunction(
    x -> ( (x[1]/R2)*(x[1]/(R2)) + (x[2]/R)*(x[2]/R) ) - 1.0,
    x -> VectorValue(2.0*(x[1]/(R2*R2)),2.0*(x[2]/(R*R))) )

  domain = (-L2,L2,0,L)
  partition = (part,div(part,2)) # If part = 10, then partition = (20,10) 
  h = (domain[2]-domain[1]) / part
  model = CartesianDiscreteModel(domain,partition) 

 # Lets write down some parameters in a txt just in case 
  pVTU="./VTU/"*simulation
  mkpath(pVTU)
  pPNG="./PNG/"*simulation
  mkpath(pPNG)
  # mkpath(pPNG*"ezrin_time/") 
  # mkpath(pPNG*"ezrin_unbound_time/") 
  mkpath(pPNG*"Rac_time/") 
  mkpath(pPNG*"Rho_time/") 
  mkpath(pPNG*"Rac_time_initial/") 
  mkpath(pPNG*"Rho_time_initial/") 
 
  # Lets copy the code in the output folder to be able to check code used for each simulation
  cp(@__FILE__, pPNG*split(@__FILE__, "/")[end],force=true)


  order = 1
  #Starting Boundary conditions
  # x₀ = 0.0
  # xₗ = 0.0
  # diri_x(p,x₀,xₗ) = p[1] < 0 ? x₀ : xₗ # x₀ on negative x coordinate, xₗ otherwise  

  v₀ = 0.0
  vₗ = 0.0
  diri_v(p,v₀,vₗ) = p[1] < 0 ? v₀ : vₗ # v₀ on negative x coordinate, vₗ otherwise  

  degree = 3 # Better for algoim quadratures
  # triangulate FE space
  Ω = Triangulation(model)
  squad = Quadrature(algoim,φ,degree,phase=CUT)
  Ωᶜ,dΓ,_ = TriangulationAndMeasure(Ω,squad) # dΓ  = quadrature rule on sphere, Ωᶜ = cut cells
  dΩᶜ = Measure(Ωᶜ,2*order)                  # dΩᶜ = quadrature rule on cut cells (for stabilisation)
  nΓ = normal(φ,Ω)

  # Creating FE space
  reffe = ReferenceFE(lagrangian,Float64,order)
  WD0 = TestFESpace(Ωᶜ,reffe;dirichlet_tags=["tag_5"])
#  X = TrialFESpace(WD0,p->diri_x(p,x₀,xₗ))
  V = TrialFESpace(WD0,p->diri_v(p,v₀,vₗ))

  #Neumann and defining trial spaces and test space for MCA_b 
  Q0 = TestFESpace(Ωᶜ,reffe)
  RHO =  TrialFESpace(Q0)
  RHO0 =  TrialFESpace(Q0)
  Rac =  TrialFESpace(Q0)
  Rho =  TrialFESpace(Q0)

  #Building the vectors used for introducing opto influence as an increase in α and β
  arclength(x) = R2 * atan(x[2],-x[1]) # Arc length for sphere
  α₀opto(x) = α₀
  β₀opto(x) = β₀

  α₀opto2(x) = α₀ + αopto * exp( -0.5 * ( arclength(x)-π*R2 )^2 / ((wrac)^2) )
  β₀opto2(x) = β₀ + βopto * exp( -0.5 * ( arclength(x))^2       / ((wrac)^2) )
  
  γ₀ = 0.1  / h # TODO: Eric reviews the scaling with h
  γ₀v = 0.1 * η / h # TODO: Eric reviews the scaling with h
  # γ₀x = 0.1 * k / h # TODO: Eric reviews the scaling with h
  γ₀M =  0.1  / h 
  γ₀R =  0.1  / h 
  m₀opto(u,v) = ∫( u*v )dΓ
  s₀opto(u,v) = ∫( γ₀*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ
  s₀v(u,v) = ∫( γ₀v*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ
  # s₀x(u,v) = ∫( γ₀x*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ
  s₀MCA(u,v) = ∫( γ₀M*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ
  s₀R(u,v) = ∫( γ₀R*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ

  A₀opto(u,v) = m₀opto(u,v) + s₀opto(u,v)
  bα₀opto(v) = m₀opto(α₀opto,v)
  bβ₀opto(v) = m₀opto(β₀opto,v)
  bα₀opto2(v) = m₀opto(α₀opto2,v)
  bβ₀opto2(v) = m₀opto(β₀opto2,v)

  op_α₀ = AffineFEOperator(A₀opto,bα₀opto,Rac,Q0)
  op_β₀ = AffineFEOperator(A₀opto,bβ₀opto,Rho,Q0)

  α₀v = solve(op_α₀)
  β₀v = solve(op_β₀)

  writevtk(Ωᶜ,"tmp",cellfields=["a"=>α₀v,"b"=>β₀v,"f"=>φ.φ])

  #defining Starting conditions for some variables and dummy variables to be able to build the equations
  uh_ezrin_b = interpolate_everywhere(kon*M0/(π*R2)/(koff+kon),RHO)
  uh_ezrin_b_old = uh_ezrin_b
  uh_ezrin_u = interpolate_everywhere(koff*M0/(π*R2)/(koff+kon),RHO0)
  uh_ezrin_u_old = uh_ezrin_u
  uh_v = zero(V)
  λ=0
  
  uh_rac = interpolate_everywhere(0.0,Rac) #0.95
  uh_rho = interpolate_everywhere(4.0,Rho) #0.62
  uh_rac_old = uh_rac
  uh_rho_old = uh_rho
  sum_uh_ezrin_u = ∑(∫(uh_ezrin_u)dΓ)
  sum_uh_ezrin_b = ∑(∫(uh_ezrin_b)dΓ)
  Minitial = sum_uh_ezrin_b + sum_uh_ezrin_b

  print("Start1 sum(MCA_b) "*string(sum_uh_ezrin_b)*", sum(MCA_u) "*string(sum_uh_ezrin_u)*", and sum(MCA_b+MCA_u) "*string(Minitial)*"\n")
  
  # Extract quadrature points and arc length array at every cell
  xΓ = dΓ.quad.cell_point.values
  xΓ = lazy_map(Reindex(xΓ),dΓ.quad.cell_point.ptrs)     # 2D array of xΓ (1 array per cell)
  alenΓ = lazy_map(Broadcasting(x->atan(x[2],-x[1])),xΓ) # Following cell order
  flat_xΓ = vcat(xΓ...)                                  # 1D "flattened" array of xΓ
  flat_alenΓ = vcat(alenΓ...)
  num_qpoints = length(flat_xΓ)

  # variables to store temporal information that we would like to plot later 
  vt = zeros(trunc(Int,T/Δt)+1,num_qpoints)  
  ract = zeros(trunc(Int,T/Δt)+1,num_qpoints)
  rhot = zeros(trunc(Int,T/Δt)+1,num_qpoints)
  ezrin_bt = zeros(trunc(Int,T/Δt)+1,num_qpoints)
  ezrin_ut = zeros(trunc(Int,T/Δt)+1,num_qpoints)
  
  #we write down the weak form of the membrane equation
  # to do backward eurler for the time evolwe give gridap a so-called Mass Term for a 
   a(k,x,w) = ∫( ( k * (∇ᵈ(x,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ
  m(MCA_b,Δt,x,w) = ∫( ( (χ*MCA_b) * (x*w) / Δt )*y )dΓ

  mMCA,aezrin_b,bezrin_b,aezrin_u,bezrin_u = MCA_bound_unbound_weak_forms(
    Δt,kon,koff,R2,D,nΓ,dΓ)


  #Now for v
  aᵥ(MCA_b,v,w) =  ∫( ( η * (∇ᵈ(v,nΓ)⋅∇ᵈ(w,nΓ)) + 0.1*(v*w) )*y )dΓ + ∫( ( (χ*MCA_b) * (v*w) )*y )dΓ
  bᵥ(w,uh_ezrin_b,uh_rho) =  ∫( ( σₐ₀*(w*(∇ᵈ(uh_rho,nΓ)⋅(TensorValue(0.0,-1.0,1.0,0.0)⋅nΓ) )) )*y )dΓ # ∫(( gradrho*w )*y )dΓ no feedback is ∫( w*∇σₐ )dΓ
  Aᵥ(v,w) = aᵥ(uh_ezrin_b,v,w) + s₀v(v,w)
  Bᵥ(w) = bᵥ(w,uh_ezrin_b,uh_rho)
  #We can now use MCA_b and x to solve v
  op_v= AffineFEOperator(Aᵥ,Bᵥ,V,WD0)

  #SOLVING X AT t=0
  # aₓ_0(x,w) = aₓ(0,x,w) 
  # b_0(w) =  bₓ(0,0,w)
  # op_x = AffineFEOperator(aₓ_0,b_0,X,WD0)
  # uh_x = solve(op_x)
  # uh_x_old =  uh_x

  #SOLVE MCA_b AT t=0 vien initial velocity zero
  Aezrin_b(MCA_b,w) = aezrin_b(MCA_b,uh_v,w) + s₀MCA(MCA_b,w)
  Bezrin_b(w) = bezrin_b(w,uh_ezrin_u,uh_ezrin_b_old,λ)
  op_ezrin_b= AffineFEOperator(Aezrin_b,Bezrin_b,RHO,Q0)
  uh_ezrin_b=solve(op_ezrin_b)
  uh_ezrin_b_old=uh_ezrin_b

  #SOLVE MCA_u AT t=0
  Aezrin_u(MCA_u,w) = aezrin_u(MCA_u,w) + s₀MCA(MCA_u,w)
  Bezrin_u(w) = bezrin_u(w,uh_ezrin_b,uh_ezrin_u_old,λ)
  op_ezrin_u = AffineFEOperator(Aezrin_u,Bezrin_u,RHO0,Q0)
  uh_ezrin_u = solve(op_ezrin_u)
  uh_ezrin_u_old = uh_ezrin_u

  #SOLVING V AT t=0 using x, xold and rho
  uh_v=zero(V)
  
  #define the production rates as vectors so that we can have heterogeneous α and β due to mechanics
  # α₀v = interpolate_everywhere(α₀,Rac) #zeros(num_qpoints) α₀v[:] .= α₀
  # β₀v = interpolate_everywhere(β₀,Rho)#zeros(num_qpoints) β₀v[:] .= β₀

  # #we define the tension for a spring
  perm=sortperm(flat_alenΓ) # Permutation to order by increasing arclength
  flat_alenΓ = R2*flat_alenΓ[perm]

  vt[1,:] = vcat(lazy_map(uh_v,xΓ)...)
  # xt[1,:] = vcat(lazy_map(uh_x,xΓ)...)
  vt[1,:] = vt[1,perm]
  # xt[1,:] = xt[1,perm]

# Computing tension
# Deleted

  # #storing information of how alpha an beta behave spatially over time
  # for j in 1:1:(partition)
  #   αt[1,j]= α₀v[j] + α*(0.5*(1-tanh((get_free_dof_values(uh_ezrin_b)[j])/rho0-ezrin_bth/rho0)))
  # end
  # βt[1,:] = β₀v .+ β*threshold(ten,sig0,tenth)

  a_rac, b_rac, a_rho, b_rho = rac_rho_weak_forms(
    Δt,dᵃ,dᵇ,Drac,Drho,α,β,α₀v,β₀v,rho0,ezrin_bth,sig0,
    tenth,a_t,b_t,uh_rac_old,uh_rho_old,nΓ,dΓ,mMCA)
 
  #SOLVE Rac AT t=0
  Arac(rac,w) = a_rac(rac,w,uh_rho,uh_ezrin_b,α₀v) + s₀R(rac,w)
  Brac(w) = b_rac(w,uh_rho,uh_ezrin_b,α₀v,uh_rac_old)
  op_rac = AffineFEOperator(Arac,Brac,Rac,Q0)
  uh_rac = solve(op_rac)
  uh_rac_old = uh_rac
  #SOLVE Rho AT t=0
  Arho(rho,w) = a_rho(rho,w,uh_rac,β₀v) + s₀R(rho,w)
  Brho(w) = b_rho(w,uh_rac,β₀v,uh_rho_old)
  op_rho = AffineFEOperator(Arho,Brho,Rho,Q0)
  uh_rho = solve(op_rho) 
  uh_rho_old = uh_rho

  dummyx0=0

  # op_α₀ = AffineFEOperator(A₀opto,bα₀opto2,Rac,Q0)
  # op_β₀ = AffineFEOperator(A₀opto,bβ₀opto2,Rho,Q0)

  # α₀v = solve(op_α₀)
  # β₀v = solve(op_β₀)
  # give steady state as initial conditions for rac and rho 
  for ti in 1:30 
    
    op_rho = AffineFEOperator(Arho,Brho,Rho,Q0)
    uh_rho = solve(op_rho)
    uh_rho_old = uh_rho
    op_rac = AffineFEOperator(Arac,Brac,Rac,Q0)
    uh_rac = solve(op_rac)
    uh_rac_old = uh_rac 
    ractt = vcat(lazy_map(uh_rac,xΓ)...)
    rhott = vcat(lazy_map(uh_rho,xΓ)...)
    ractt[:] = ractt[perm]
    rhott[:] = rhott[perm]
    plotting("rac",ractt[:],pPNG*"Rac_time_initial/","$ti")
    plotting("rho",rhott[:],pPNG*"Rho_time_initial/","$ti")
  end
  # threshold value to start protrusion, 1.3 times the initial condition 
  racaux=vcat(lazy_map(uh_rac,xΓ)...)
  racaux=racaux[perm]
  racth = 1.3*racaux[end] 
  
  ract[1,:] = vcat(lazy_map(uh_rac,xΓ)...)
  rhot[1,:] = vcat(lazy_map(uh_rho,xΓ)...)
  ract[1,:] = ract[1,perm]
  rhot[1,:] = rhot[1,perm]
  
  i = 0
  t=0
  writevtk(Ωᶜ,pVTU*"VTU$i",cellfields=["v"=>uh_v,"ezrin_b"=>uh_ezrin_b,"ezrin_u"=>uh_ezrin_u,"rac"=>uh_rac,"rho"=>uh_rho,"f"=>φ.φ]) 
 for ti in t₀:Δt:(T-Δt)
    #HERE WE DEFINE WHETHER THE CODE IS FRONT TO BACK OR BACK TO FRONT, DEPENDING IN WHERE WE ACTIVATE OPTO
    if t==topto 
      op_α₀ = AffineFEOperator(A₀opto,bα₀opto2,Rac,Q0)
      op_β₀ = AffineFEOperator(A₀opto,bβ₀opto2,Rho,Q0)

      α₀v = solve(op_α₀)
      β₀v = solve(op_β₀)
      α₀vaux = vcat(lazy_map(α₀v,xΓ)...)
      α₀vaux[:] = α₀vaux[perm]
    
    end
    if t==3*topto #at time=3*topto a while the input dies down
      op_α₀ = AffineFEOperator(A₀opto,bα₀opto,Rac,Q0)
      op_β₀ = AffineFEOperator(A₀opto,bβ₀opto,Rho,Q0)

      α₀v = solve(op_α₀)
      β₀v = solve(op_β₀)
    end
    i1 = ∑(∫(uh_ezrin_b)dΓ)
    i2 = ∑(∫(uh_ezrin_u)dΓ)
    λ = conservation(i2,i1,Minitial)
    i3 = i1+i2
    i4 = trunc(t)
    i = i + 1
    t = t + Δt  

    #boundary conditions for v and x
    racaux=vcat(lazy_map(uh_rac,xΓ)...)
    racaux=racaux[perm] 
    rac1=racaux[end] 
   # ten1=tenaux[end]
    #vCTErac = vCTE*threshold2(rac1,rac0,1.3*racaux[1])#velocity polimerization
    #vₗ = 0.5
    #we introduce a slight relaxation for the membrane, decreases 2% x at the Boundary condition only
    dummyx0=dummyx0*0.9+Δt*vₗ
    #xₗ = dummyx0
    #X = TrialFESpace(WD0,p->diri_x(p,x₀,xₗ))

    sum_uh_a = ∑(∫(uh_rac)dΓ)
    sum_uh_b = ∑(∫(uh_rho)dΓ)

    @info "Time step $i/$nΔt, time $i4, sum(MCA_b+MCA_u) $i3 sum Rac $sum_uh_a sum rho $sum_uh_b"
    #Mt[i]=i3

    ezrin_bt[i+1,:] = vcat(lazy_map(uh_ezrin_b,xΓ)...)
    ezrin_ut[i+1,:] = vcat(lazy_map(uh_ezrin_u,xΓ)...)
    ract[i+1,:] = vcat(lazy_map(uh_rac,xΓ)...)
    rhot[i+1,:] = vcat(lazy_map(uh_rho,xΓ)...)
    ezrin_bt[i+1,:] = ezrin_bt[i+1,perm]
    ezrin_ut[i+1,:] = ezrin_ut[i+1,perm]
    ract[i+1,:] = ract[i+1,perm]
    rhot[i+1,:] = rhot[i+1,perm]

    # Updating v to solve ezrin_b and x
    # A(x,w) = m(uh_ezrin_b,Δt,x,w) + a(k,x,w) + s₀x(x,w)
    # B(w) = m(uh_ezrin_b,Δt,uh_x,w) + bₓ(uh_ezrin_b,uh_v,w) 
    # op_x = AffineFEOperator(A,B,X,WD0)
    # uh_x = solve(op_x)

    # op_ezrin_b= AffineFEOperator(Aezrin_b,Bezrin_b,RHO,Q0)
    # uh_ezrin_b=solve(op_ezrin_b)
    # uh_ezrin_b_old=uh_ezrin_b
  
    # op_ezrin_u= AffineFEOperator(Aezrin_u,Bezrin_u,RHO0,Q0)
    # uh_ezrin_u=solve(op_ezrin_u)
    # uh_ezrin_u_old=uh_ezrin_u


    V = TrialFESpace( WD0,p->diri_v(p, v₀ , vₗ) )
    op_v = AffineFEOperator(Aᵥ,Bᵥ,V,WD0)
    uh_v = solve(op_v)# uh_v = interpolate_everywhere(f3,WD0) #
    #f4(x) = vCTE - vCTE/(R*π)*arclength(x)
    #uh_v = interpolate_everywhere(f4,WD0) 

    #if t<10*topto 
    # a_sum = sum(pa) 
  #  Arac(rac,w) = a_rac(rac,w,uh_rho,uh_ezrin_b,α₀v)
  #  Brac(w) = b_rac(w,uh_rho,uh_ezrin_b,α₀v,uh_rac_old)
      op_rac = AffineFEOperator(Arac,Brac,Rac,Q0)
      uh_rac = solve(op_rac)
      uh_rac_old = uh_rac
      
      # b_sum = sum(paa) 
      #Arho(rho,w) = a_rho(rho,w,ten,uh_rac,β₀v)
      #Brho(w) = b_rho(w,ten,uh_rac,β₀v,uh_rho_old)
      op_rho = AffineFEOperator(Arho,Brho,Rho,Q0)
      uh_rho = solve(op_rho)
      uh_rho_old = uh_rho
    #end

    #updating x_old per time derivarive
    # uh_x_old =  uh_x
    writevtk(Ωᶜ,pVTU*"VTU$i",cellfields=["v"=>uh_v,"ezrin_b"=>uh_ezrin_b,"ezrin_u"=>uh_ezrin_u,"rac"=>uh_rac,"rho"=>uh_rho,"f"=>φ.φ]) 
 
    #plotting("MCA_b",ezrin_bt[i+1,:],pPNG*"ezrin_time/","$i")
    plotting("rac",ract[i+1,:],pPNG*"Rac_time/","$i")
    plotting("rho",rhot[i+1,:],pPNG*"Rho_time/","$i")
   # plotting("MCA_u",ezrin_ut[i+1,:],pPNG*"ezrin_unbound_time/","$i")



    # op_ten = AffineFEOperator(Aten,bten,Rac,Q0) 

    # ten = solve(op_ten)


    vt[i+1,:] = vcat(lazy_map(uh_v,xΓ)...)
    # xt[i+1,:] = vcat(lazy_map(uh_x,xΓ)...)
    vt[i+1,:] = vt[i+1,perm]
    # xt[i+1,:] = xt[i+1,perm]
    # tenaux = vcat(lazy_map(ten,xΓ)...)
    # tensiont[i+1,:] = tenaux[perm]

    
  end 
  plots_run_singlet(nΔt,h,vt,ezrin_bt,ezrin_ut,ract,rhot,pPNG,α,
   β, dᵃ, dᵇ,αopto,βopto,topto,num_qpoints,L,Δt,T,
    η, σₐ₀ , ten0, rho0,α₀,β₀,flat_alenΓ,part) 
  print("finish line")
  return ezrin_bt,vt,ract,rhot
end