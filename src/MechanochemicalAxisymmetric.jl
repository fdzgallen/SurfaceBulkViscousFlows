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

function conservation(sMCAu,sMCAb,Minitial)
  return 0.1*(Minitial-(sMCAu+sMCAb))
end
 
function threshold(x,x₀,xth)
  return  (0.5 * (tanh.(x/x₀ .- xth/x₀).+1)) 
end

#function to run a single simulation with a few given parameters
function run_mechanochemical_axisymmetric(χ,λ⁻²,η,T,Δt,part,
    L,simulation,wrac,αopto,βopto,kon,koff,M0,α₀,β₀,k,D,
    σₐ₀,λᵇ,Drac,Drho,mTorc,rho0,ten0,a_t,b_t,α,β,dᵃ,dᵇ,sig0,tenth,λʳᴬ)

  # Time discretisation parameters
  t₀  = 0.0
  t   = t₀  
  nΔt = trunc(Int,T/Δt)

  # Level set function implicitly describing a unit sphere
  R = 1.0
  φ = AlgoimCallLevelSetFunction(
    x -> ( (x[1]/R)*(x[1]/R) + (x[2]/R)*(x[2]/R) ) - 1.0,
    x -> VectorValue(2.0*(x[1]/(R*R)),2.0*(x[2]/(R*R))) )

  domain = (-L,L,0,L)
  partition = (part,div(part,2)) # If part = 10, then partition = (20,10) 
  h = (domain[2]-domain[1]) / part
  model = CartesianDiscreteModel(domain,partition) 

 # Lets write down some parameters in a txt just in case 
  pVTU="./VTU/"*simulation
  mkpath(pVTU)
  pPNG="./PNG/"*simulation
  mkpath(pPNG)
  mkpath(pPNG*"ezrin_time/") 
  mkpath(pPNG*"ezrin_unbound_time/") 
  mkpath(pPNG*"Rac_time/") 
  mkpath(pPNG*"Rho_time/") 
 
  # Lets copy the code in the output folder to be able to check code used for each simulation
  cp(@__FILE__, pPNG*split(@__FILE__, "/")[end],force=true)


  order = 1
  #Starting Boundary conditions
  x₀ = 0
  xₗ = 0
  diri_x(p) = p[1] < 0 ? x₀ : xₗ # x₀ on negative x coordinate, xₗ otherwise  

  v₀ = 0
  vₗ = 0
  diri_v(p) = p[1] < 0 ? v₀ : vₗ # v₀ on negative x coordinate, vₗ otherwise  

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
  X = TrialFESpace(WD0,diri_x)
  V = TrialFESpace(WD0,diri_v)

  #Neumann and defining trial spaces and test space for MCA_b 
  Q0 = TestFESpace(Ωᶜ,reffe)
  RHO =  TrialFESpace(Q0)
  RHO0 =  TrialFESpace(Q0)
  Rac =  TrialFESpace(Q0)
  Rho =  TrialFESpace(Q0)

  #Building the vectors used for introducing opto influence as an increase in α and β
  arclength(x) = R * atan(x[2],-x[1]) # Arc length for sphere
  α₀opto(x) = αopto * exp( -0.5 *   arclength(x)^2       / ((wrac)^2))
  β₀opto(x) = βopto * exp( -0.5 * ( arclength(x)-π*R )^2 / ((wrac)^2))
  
  γ₀ = 0.1 / h # TO-DO: Eric REVIEW's the scaling with h
  m₀opto(u,v) = ∫( u*v )dΓ
  s₀opto(u,v) = ∫( γ₀*((nΓ⋅∇(u))⊙(nΓ⋅∇(v))) )dΩᶜ

  A₀opto(u,v) = m₀opto(u,v) + s₀opto(u,v)
  bα₀opto(v) = m₀opto(α₀opto,v)
  bβ₀opto(v) = m₀opto(β₀opto,v)

  op_α₀ = AffineFEOperator(A₀opto,bα₀opto,Rac,Q0)
  op_β₀ = AffineFEOperator(A₀opto,bβ₀opto,Rho,Q0)

  α₀v = solve(op_α₀)
  β₀v = solve(op_β₀)

  writevtk(Ωᶜ,"tmp",cellfields=["a"=>α₀v,"b"=>β₀v,"f"=>φ.φ])

  #defining Starting conditions for some variables and dummy variables to be able to build the equations
  uh_MCAb = interpolate_everywhere(kon*M0/part/(koff+kon),RHO)
  uh_MCAb_old = uh_MCAb
  uh_MCAu = interpolate_everywhere(koff*M0/part/(koff+kon),RHO0)
  uh_MCAu_old = uh_MCAu
  uh_x = zero(X)
  uh_x_old = uh_x
  uh_v = zero(V)
  f1(x) = 0.01*arclength(x)
  f2(x) = α₀+0.5*exp(-arclength(x)^2/wrac^2)
  uh_rac = interpolate_everywhere(0.0,Rac) #0.95
  uh_rho = interpolate_everywhere(0.62,Rho) #0.62
  uh_rac_old = uh_rac
  uh_rho_old = uh_rho
  sum_uh_MCAu = ∑(∫(uh_MCAu)dΓ)
  sum_uh_MCAb = ∑(∫(uh_MCAb)dΓ)
  Minitial = sum_uh_MCAu + sum_uh_MCAb
  
  print("Start1 sum(MCA_b) "*string(sum_uh_MCAb)*", sum(MCA_u) "*string(sum_uh_MCAu)*", and sum(MCA_b+MCA_u) "*string(Minitial)*"\n")
  
  # Extract quadrature points and arc length array at every cell
  xΓ = dΓ.quad.cell_point.values
  xΓ = lazy_map(Reindex(xΓ),dΓ.quad.cell_point.ptrs)     # 2D array of xΓ (1 array per cell)
  alenΓ = lazy_map(Broadcasting(x->atan(x[2],-x[1])),xΓ) # Following cell order
  flat_xΓ = vcat(xΓ...)                                  # 1D "flattened" array of xΓ
  flat_alenΓ = vcat(alenΓ...)
  num_qpoints = length(flat_xΓ)

  # variables to store temporal information that we would like to plot later
  tensiont = zeros(trunc(Int,T/Δt)+1,num_qpoints-1)
  vt = zeros(trunc(Int,T/Δt)+1,num_qpoints)
  xt = zeros(trunc(Int,T/Δt)+1,num_qpoints)
  αt = zeros(trunc(Int,T/Δt)+1,num_qpoints)
  βt = zeros(trunc(Int,T/Δt)+1,num_qpoints)
  Mt = zeros(trunc(Int,T/Δt))
  ract = zeros(trunc(Int,T/Δt)+1,num_qpoints)
  rhot = zeros(trunc(Int,T/Δt)+1,num_qpoints)
  #We start wriiting down the values of some parameters at time=0
  λ = 0.0
  MCAbt = zeros(trunc(Int,T/Δt)+1,num_qpoints)
  MCAbt[1,:] = vcat(lazy_map(uh_MCAb,xΓ)...)
  MCAut = zeros(trunc(Int,T/Δt)+1,num_qpoints)
  MCAut[1,:] = vcat(lazy_map(uh_MCAu,xΓ)...)  
  
  #we write down the weak form of the membrane equation
  # to do backward eurler for the time evolwe give gridap a so-called Mass Term for a 
  a(k,x,w) = ∫( ( k * (∇ᵈ(x,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ
  m(MCA_b,Δt,x,w) = ∫( ( (χ*MCA_b) * (x*w) / Δt )*y )dΓ
  
  aₓ(MCA_b,x,w) = ∫( ( (χ*MCA_b) * (x*w) / Δt )*y )dΓ + ∫( ( k * (∇ᵈ(x,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ
  bₓ(MCA_b,v,w) = ∫( ( χ*(MCA_b*v)*w )*y )dΓ
  
  # Now for MAC bound
  # mass term for the temporal evolution MCA_b
  mMCA(Δt,MCA_b,w) = ∫( ( (MCA_b*w)/Δt )*y )dΓ
  aMCAb(MCA_b,v,w) = mMCA(Δt,MCA_b,w) + ∫( ( 0.001 * (∇ᵈ(MCA_b,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ + 
    ∫( ( w * ( v * (∇ᵈ(MCA_b,nΓ)⋅VectorValue(1.0,1.0)) + MCA_b * (∇ᵈ(v,nΓ)⋅VectorValue(1.0,1.0)) ) )*y )dΓ + 
    ∫( ( koff * ( MCA_b * w )*y ) )dΓ + ∫( ( λᵇ * ( ( MCA_b*MCA_b*MCA_b ) * w ) )*y )dΓ
  bMCAb(w,MCA_u) = ∫( ( kon * ( MCA_u * w ) + λ/(π*R) * ( kon / (kon+koff) ) )*y )dΓ + mMCA(Δt,uh_MCAb_old,w)  

  #Now for MCA unbound
  # mass term for the temporal evolution MCA_b0
  aMCAu(MCA_u,x,x_old,w) = mMCA(Δt,MCA_u,w) + 
    ∫( ( kon * ( MCA_u * w ) )*y )dΓ + ∫( ( D * ( ∇ᵈ(MCA_u,nΓ)⋅∇ᵈ(w,nΓ) ) )*y )dΓ + 
    ∫( ( w * ( ( (x-x_old) / Δt ) * (∇ᵈ(MCA_u,nΓ)⋅VectorValue(1.0,1.0)) + 
          MCA_u * ( (∇ᵈ(x,nΓ)⋅VectorValue(1.0,1.0)) - (∇ᵈ(x_old,nΓ)⋅VectorValue(1.0,1.0)) ) / Δt ) )*y )dΓ  
  bMCAu(w,MCA_b) = ∫( ( koff * (MCA_b*w) + λ/(π*R) * (koff/(kon+koff)) )*y )dΓ + mMCA(Δt,uh_MCAu_old,w) 

  #Now for v
  aᵥ(MCA_b,v,w) = ∫( ( η * (∇ᵈ(v,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ + ∫( ( (χ*MCA_b) * (v*w) )*y )dΓ
  bᵥ(w,uh_rho) = m(uh_MCAb,Δt,uh_x,w) - m(uh_MCAb,Δt,uh_x_old,w) + 
    ∫( ( σₐ₀*(w*(∇ᵈ(uh_rho,nΓ)⋅VectorValue(1.0,1.0))) )*y )dΓ #no feedback is ∫( w*∇σₐ )dΓ
  Aᵥ(v,w) = aᵥ(uh_MCAb,v,w)
  Bᵥ(w) = bᵥ(w,uh_rho)
  #We can now use MCA_b and x to solve v
  op_v= AffineFEOperator(Aᵥ,Bᵥ,V,WD0)

  #SOLVING X AT t=0
  aₓ_0(x,w) = aₓ(0,x,w)
  b_0(w) =  bₓ(0,0,w)
  op_x = AffineFEOperator(aₓ_0,b_0,X,WD0)
  uh_x = solve(op_x)
  uh_x_old =  uh_x

  #SOLVE MCA_b AT t=0 vien initial velocity zero
  AMCAb(MCA_b,w) = aMCAb(MCA_b,uh_v,w) 
  BMCAb(w) = bMCAb(w,uh_MCAu)
  op_MCAb= AffineFEOperator(AMCAb,BMCAb,RHO,Q0)
  uh_MCAb=solve(op_MCAb)
  uh_MCAb_old=uh_MCAb

  #SOLVE MCA_u AT t=0
  AMCAu(MCA_u,w) = aMCAu(MCA_u,uh_x,uh_x_old,w)
  BMCAu(w) = bMCAu(w,uh_MCAb_old)
  op_MCAu = AffineFEOperator(AMCAu,BMCAu,RHO0,Q0)
  uh_MCAu = solve(op_MCAu)
  uh_MCAu_old = uh_MCAu

  #SOLVING V AT t=0 using x, xold and rho
  uh_v=zero(V)
  
  #define the production rates as vectors so that we can have heterogeneous α and β due to mechanics
  α₀v = zeros(num_qpoints)
  α₀v[:] .= α₀
  β₀v = zeros(num_qpoints)
  β₀v[:] .= β₀

  #we define the tension for a spring
  perm=sortperm(flat_alenΓ) # Permutation to order by increasing arclength
  lx=vcat(lazy_map(uh_x,xΓ)...)
  lx=lx[perm]
  x1=circshift(lx,1)
  tension=k*(lx-x1)/h 
  splice!(tension,  1) 

  vt[1,:] = vcat(lazy_map(uh_v,xΓ)...)
  xt[1,:] = vcat(lazy_map(uh_x,xΓ)...)    
  tensiont[1,:] = tension[:]
  #ten will be used in the computation of the model, we instill the Boundary Conditions and remove negative values
  ten = zeros(num_qpoints, 1)
  ten[1:end-1] = tension[:]
  ten[end-1] = tension[end]
  ten[end] = tension[end]
  
  # #storing information of how alpha an beta behave spatially over time
  # for j in 1:1:(partition)
  #   αt[1,j]= α₀v[j] + α*(0.5*(1-tanh((get_free_dof_values(uh_MCAb)[j])/rho0-MCAbth/rho0))) 
  # end
  # βt[1,:] = β₀v .+ β*threshold(ten,sig0,tenth) 

  #DEFINING the equations for Rac and Rho
  a_rac(rac,w,rho,MCA_b) = (1/dᵃ) * mMCA(Δt,rac,w) + 
    ∫( ( Drac * (∇ᵈ(rac,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ + 
    ∫( ( w*rac )*y )dΓ + ∫( ( w*rac * ( α₀v/(1+rho*rho) + α*(0.5*(1-tanh∘(MCA_b/rho0-MCAbth/rho0))) / (1+rho*rho) ) )*y )dΓ + 
    ∫( ( mTorc*w*rac*threshold(ten,ten0,10) )*y )dΓ 
  b_rac(w,rho,MCA_b) = (1/dᵃ) * mMCA(Δt,uh_rac_old,w) + 
    ∫( ( w*(a_t)*(α₀v/(1+rho*rho) + α*(0.5*(1-tanh∘(MCA_b/rho0-MCAbth/rho0)))/(1+rho*rho)) )*y )dΓ  
 
  a_rho(rho,w,ten,rac) = (1/dᵇ)*mMCA(Δt,rho,w) + 
    ∫( ( Drho * (∇ᵈ(rho,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ + 
    ∫( ( w*rho )*y )dΓ + 
    ∫( ( w*rho * ( β₀v/(1+rac*rac) + β*threshold(ten,sig0,tenth)/(1+rac*rac) ) )*y )dΓ +
    ∫( ( λʳᴬ*((rho*rho*rho)*w) )*y )dΓ  
  b_rho(w,ten,rac) = (1/dᵇ)*mMCA(Δt,uh_rho_old,w) + ∫( ( w*b_t*(β₀v/(1+rac*rac) + β*threshold(ten,sig0,tenth)/(1+rac*rac)) )*y )dΓ  
 
  #SOLVE Rac AT t=0
  Arac(rac,w) = a_rac(rac,w,uh_rho,uh_MCAb)
  Brac(w) = b_rac(w,uh_rho,uh_MCAb)
  op_rac = AffineFEOperator(Arac,Brac,Rac,Q0)
  uh_rac = solve(op_rac)
  uh_rac_old = uh_rac
  #SOLVE Rho AT t=0
  Arho(rho,w) = a_rho(rho,w,ten,uh_rac)
  Brho(w) = b_rho(w,ten,uh_rac)
  op_rho = AffineFEOperator(Arho,Brho,Rho,Q0)
  uh_rho = solve(op_rho)
  uh_rho_old = uh_rho

  i = 0
  t=0
  dummyx0=0
 #give steady state as initial conditions for rac and rho 
  for ti in 1:250
    op_rac = AffineFEOperator(Arac,Brac,Rac,Q0)
    uh_rac = solve(op_rac)
    uh_rac_old = uh_rac
    op_rho = AffineFEOperator(Arho,Brho,Rho,Q0)
    uh_rho = solve(op_rho)
    uh_rho_old = uh_rho
  end
  #threshold value to start protrusion, 1.3 times the initial condition
  racth = 1.3*get_free_dof_values(uh_rac)[1] 

  
  ract[1,:] = vcat(lazy_map(uh_rac,xΓ)...)
  rhot[1,:] = vcat(lazy_map(uh_rho,xΓ)...)
  
  writevtk(Ω,pVTU*"VTU$i",cellfields=["x"=>uh_x,"v"=>uh_v,"rho"=>uh_MCAb,"rho0"=>uh_MCAu]) 
  for ti in t₀:Δt:(T-Δt)
    #HERE WE DEFINE WHETHER THE CODE IS FRONT TO BACK OR BACK TO FRONT, DEPENDING IN WHERE WE ACTIVATE OPTO
#     if t==topto
#       α₀v[:] .= α₀ .+ α₀opto[:]
#       β₀v[:] .= β₀ .+ β₀opto[:]
#     end
#     if t==3*topto #at time=3*topto a while the input dies down
#       α₀v[:] .= α₀
#       β₀v[:] .= β₀ 
#     end
#     i1 = sum(∫(uh_MCAb)dΓ)
#     i2 = sum(∫(uh_MCAu)dΓ)
#     λ = conservation(i2,i1,Minitial)
#     i3 = i1+i2
#     i4 = trunc(t)
#     i = i + 1
#     t = t + Δt  

#     #boundary conditions for v and x
#     rac1=get_free_dof_values(uh_rac)[1]
#     vCTErac = vCTE*threshold(rac1,rac0,racth)#velocity polimerization
#     v₀ = vCTErac/(1+ten[1]*ten[1]/ten0)
#     #we introduce a slight relaxation for the membrane, decreases 2% x at the Boundary condition only
#     dummyx0=dummyx0*0.98+Δt*vCTErac/(1+ten[1]*ten[1]/ten0)
#     x₀ = dummyx0
#     X = TrialFESpace(WD0,[x₀,xₗ])

#     @info "Time step $i/$nΔt, time $i4, sum(MCA_b+MCA_u) $i3"
#     Mt[i]=i3

#     MCAbt[i+1,:] = get_free_dof_values(uh_MCAb)[:]
#     MCAut[i+1,:] = get_free_dof_values(uh_MCAu)[:]
#     ract[i+1,:] = get_free_dof_values(uh_rac)[:]
#     rhot[i+1,:] = get_free_dof_values(uh_rho)[:]

#     # Updating v to solve MCAb and x
#     A(x,w) = m(uh_MCAb,Δt,x,w) + a(k,x,w)
#     B(w) = m(uh_MCAb,Δt,uh_x,w) + bₓ(uh_MCAb,uh_v,w) 
#     op_x = AffineFEOperator(A,B,X,WD0)
#     uh_x = solve(op_x)

#     po=get_free_dof_values(uh_MCAb)
#     pa=get_free_dof_values(uh_rac)
#     paa=get_free_dof_values(uh_rho)
#     poo = get_free_dof_values(uh_MCAu)

#     op_MCAb= AffineFEOperator(AMCAb,BMCAb,RHO,Q0)
#     uh_MCAb=solve(op_MCAb)
#     uh_MCAb_old=uh_MCAb
  
#     op_MCAu= AffineFEOperator(AMCAu,BMCAu,RHO0,Q0)
#     uh_MCAu=solve(op_MCAu)
#     uh_MCAu_old=uh_MCAu
 
#     V = TrialFESpace(WD0,[v₀,vₗ])
#     op_v= AffineFEOperator(Aᵥ,Bᵥ,V,WD0)
#     uh_v=solve(op_v)

#     a_sum = sum(pa) 
#     op_rac = AffineFEOperator(Arac,Brac,Rac,Q0)
#     uh_rac = solve(op_rac)
#     uh_rac_old = uh_rac
    
#     b_sum = sum(paa) 
#     op_rho = AffineFEOperator(Arho,Brho,Rho,Q0)
#     uh_rho = solve(op_rho)
#     uh_rho_old = uh_rho

#     #updating x_old per time derivarive
#     uh_x_old =  uh_x
#     writevtk(Ω,pVTU*"VTU$i",cellfields=["x"=>uh_x,"v"=>uh_v,"rho"=>uh_MCAb,"rho0"=>uh_MCAu])

#     plotting("MCA_b",po,pPNG*"ezrin_time/","$i")
#     plotting("rac",pa,pPNG*"Rac_time/","$i")
#     plotting("rho",paa,pPNG*"Rho_time/","$i")
#     plotting("MCA_u",poo,pPNG*"ezrin_unbound_time/","$i")

#     lx=get_free_dof_values(uh_x)
#     x1=circshift(lx,1)
#     tension=k*(lx-x1)/h
#     splice!(tension,  1)
#     ten[1:end-2] = tension[:]
#     ten[end-1] = tension[end]
#     ten[end] = tension[end] 
#     vt[i+1,:] = get_free_dof_values(uh_v)[:]
#     xt[i+1,:] = get_free_dof_values(uh_x)[:]
#     tensiont[i+1,:] = tension[:]
#     βt[i+1,:] = β₀v .+ β*threshold(ten,sig0,tenth) 
#     for j in 1:1:(partition)
#       αt[i+1,j] = α₀v[j] + α*(0.5*(1-tanh((get_free_dof_values(uh_MCAb)[j])/rho0-MCAbth/rho0))) 
#     end
  end 
  # plots_run(nΔt,h,vt,xt,tensiont,MCAbt,MCAut,ract,rhot,Mt,λ⁻²,pPNG,α, β, dᵃ, dᵇ,αt,βt,αopto,βopto,topto)
  # return tensiont,MCAbt,vt,ract,rhot,αt,βt
end