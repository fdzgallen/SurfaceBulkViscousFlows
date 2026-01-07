iy(x) = VectorValue( 0.0, 1.0 / x[2] ); y(x) = x[2]


function rac_rho_weak_forms2(Δt,dᵃ,dᵇ,Drac,Drho,nΓ,dΓ)
  #DEFINING the equations for Rac and Rho
  m2(Δt,A,w) = ∫( ( (A*w)/Δt )*y )dΓ
  advection(rac,w,v) =  ∫( (w*( ∇ᵈ(rac,nΓ)⋅v + rac*divᶜ(v,nΓ))  )*y )dΓ

  a_rac(rac,w,v) = (1/dᵃ) * m2(Δt,rac,w)   + ∫( ( Drac * (∇ᵈ(rac,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ + ∫( ( w*rac )*y )dΓ
  b_rac(w,rho,α₀v,rac_old) = (1/dᵃ) * m2(Δt,rac_old,w) +   ∫( ( w*(α₀v)/(1+rho*rho) )*y )dΓ   
 
  a_rho(rho,w,v) = (1/dᵇ)*m2(Δt,rho,w) + advection(rho,w,v) +
    ∫( ( Drho * (∇ᵈ(rho,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ   + ∫( ( w*rho )*y )dΓ # + ∫( ( λʳᴬ*((rho*rho*rho)*w) )*y )dΓ  
  b_rho(w,rac,β₀v,rho_old) = (1/dᵇ)*m2(Δt,rho_old,w) + ∫( ( w*(β₀v)/(1+rac*rac) )*y )dΓ 

  a_rac, b_rac, a_rho, b_rho
end



function rac_rho_weak_forms2(Δt,dᵃ,dᵇ,Drac,Drho,nΓ,dΓ,α,β)
  function threshold(x,x₀,xth)
    return  (0.5 * (tanh∘(x/x₀ - xth/x₀)+1))
  end
  #(    α,β,,,,, tenth)
  #DEFINING the equations for Rac and Rho
  m2(Δt,A,w) = ∫( ( (A*w)/Δt )*y )dΓ
  advection(rac,w,v) =  ∫( (w*( ∇ᵈ(rac,nΓ)⋅v + rac*divᶜ(v,nΓ))  )*y )dΓ

  a_rac(rac,w,v) = (1/dᵃ) * m2(Δt,rac,w)   + ∫( ( Drac * (∇ᵈ(rac,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ + ∫( ( w*rac )*y )dΓ
  b_rac(w,rho,α₀v,rac_old, MCA_b,MCAbth,rho0) = (1/dᵃ) * m2(Δt,rac_old,w) +  ∫( ( w*(α₀v + α*(1.0 - threshold(MCA_b,rho0,MCAbth)))/(1+rho*rho) )*y )dΓ    
 
  a_rho(rho,w,v) = (1/dᵇ)*m2(Δt,rho,w) +  
    ∫( ( Drho * (∇ᵈ(rho,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ   + ∫( ( w*rho )*y )dΓ # + ∫( ( λʳᴬ*((rho*rho*rho)*w) )*y )dΓ  
  b_rho(w,rac,β₀v,rho_old,ten,sig0,tenth) = (1/dᵇ)*m2(Δt,rho_old,w) + ∫( ( w*(β₀v + β*threshold(ten,sig0,tenth))/(1+rac*rac) )*y )dΓ 

  a_rac, b_rac, a_rho, b_rho
end

function rac_rho_weak_forms_conserved(Δt,dᵃ,dᵇ,Drac,Drho,nΓ,dΓ)
  #DEFINING the equations for Rac and Rho
  #TODO add y's, convert into axisymmetric
  m2(Δt,A,w) = ∫( ( (A*w)/Δt )*y )dΓ
  advection(rac,w,v) =  ∫( (w*( ∇ᵈ(rac,nΓ)⋅v + rac*divᶜ(v,nΓ))  )*y )dΓ

  a_rac(rac,w,v) = (1/dᵃ)*m2(Δt,rac,w)  + ∫( ( Drac * (∇ᵈ(rac,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ + ∫((w*rac)*y)dΓ# + advection(rac,w,v)
  b_rac(w,rho,α₀v,rac_i,rac_old) =  (1/dᵃ)*m2(Δt,rac_old,w) + ∫( w*(rac_i)*(α₀v/(1+rho*rho))*y )dΓ  
  a_rac_i(rac_i,w) = ∫(((w*(rac_i))*y))dΓ  
  b_rac_i(w,a_t,a_sum) = ∫((w*(a_t - a_sum))*y)dΓ

  a_rho(rho,w,v) = (1/dᵇ)*m2(Δt,rho,w)  + ∫( ( Drho * (∇ᵈ(rho,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ  + advection(rho,w,v) + 
  ∫((w*rho)*y)dΓ  #  +∫(λʳᴬ*((rho*rho*rho)*w))dΓ  
  b_rho(w,rac,β₀v,rho_i,rho_old) =  (1/dᵇ)*m2(Δt,rho_old,w) + ∫( (w*(rho_i)*(β₀v/(1+rac*rac)))*y)dΓ   
  a_rho_i(rho_i,w) = ∫((w*(rho_i))*y)dΓ 
  b_rho_i(w,b_t,b_sum) = ∫((w*(b_t - b_sum))*y)dΓ

  a_rac, b_rac, a_rho, b_rho, a_rac_i, b_rac_i, a_rho_i, b_rho_i
end


function cortical_flow_problem_axisymmetricOG(
    eₕ,dΩᶜ,dΓ,nΓ,γ::Float64,Pe::Float64,χ::Float64,ξ₀::Function)

  # Viscous term
  aʷ(u,v) = 
    ∫( ( εᶜ(u,nΓ)⊙εᵈ(v,nΓ) + divᶜ(u,nΓ)⋅divᶜ(v,nΓ) + 
         2*(u⋅iy)*(v⋅iy) + divᶜ(u,nΓ)*(v⋅iy) + 
         divᶜ(v,nΓ)*(u⋅iy) )*y )dΓ

  # Friction term
  aᶠ(u,v) = ∫( χ*(u⋅v)*y )dΓ

  # Activity function (relates myosin concentration to active stress)
  ξ(e) = 2.0 * e*e / ( 1.0 + e*e )

  # Active force term
  f(μ,e) = ∫( Pe * ( -(divᶜ(μ,nΓ)+μ⋅iy)*(ξ∘(e)) ) * ξ₀ )dΓ

  # Stabilisation term for velocity
  sᵘ(υ,μ) = ∫( γ * ((nΓ⋅ε(υ))⊙(nΓ⋅ε(μ))) )dΩᶜ

  # Rigid body motion and volum constraint
  RB¹ = VectorValue(1.0,0.0)
  r¹(u,ℓ) = ∫( ( u⋅(ℓ*RB¹ ) )*y )dΓ#∫( ( RB¹⋅(ℓ*u) )*y )dΓ#
  r²(u,ℓ) = ∫( ( u⋅(ℓ*nΓ ) )*y )dΓ # nΓ⋅(ℓ*u) )*y )dΓ # u⋅(ℓ*nΓ ) )*y )dΓ

  aᵛ((υ,l¹,l²),(μ,ℓ¹,ℓ²)) =
    aʷ(υ,μ) + aᶠ(υ,μ) + sᵘ(υ,μ) + 
    r¹(υ,ℓ¹) + r¹(μ,l¹) + r²(υ,ℓ²) + r²(μ,l²)
  bᵛ((μ,ℓ¹,ℓ²)) = f(μ,eₕ)

  aᵛ, bᵛ
end

function cortical_flow_problem_axisymmetric_turnover(
    ρₕ,R,eₕ,dΩᶜ,dΓ,nΓ,γ::Float64,Pe::Float64,χᵣ::Float64,χ₀::Float64,ξ₀,σₐ⁰,
  sigmaₐ⁰,  sigmaρ⁰, sigmaR⁰)

  # Viscous term
  aʷ(u,v,e) = 
    ∫( ( e*εᶜ(u,nΓ)⊙εᵈ(v,nΓ) + e*divᶜ(u,nΓ)⋅divᶜ(v,nΓ) + 
         2*e*(u⋅iy)*(v⋅iy) + e*divᶜ(u,nΓ)*(v⋅iy) + 
         e*divᶜ(v,nΓ)*(u⋅iy) )*y )dΓ

  χ(R) = (χ₀+χᵣ*R)

  # Friction term
  aᶠ(u,v,R) = ∫(χ(R)*(u⋅v)*y )dΓ

  # Activity function (relates Rho concentration to active stress) 
  function sigmaₐ(ρ,R)
      sigmaₐ = sigmaₐ⁰ + sigmaρ⁰ * ρ - sigmaR⁰ * R
      sigmaₐ > 0 ? sigmaₐ : zero(typeof(sigmaₐ))
  end
  # Active force term
  f(μ, ρ,R,e) = ∫( e*( -(divᶜ(μ,nΓ)+μ⋅iy)*(sigmaₐ∘(ρ,R))*σₐ⁰ ) * ξ₀ )dΓ 

  # Stabilisation term for velocity
  sᵘ(υ,μ) = ∫( γ * ((nΓ⋅ε(υ))⊙(nΓ⋅ε(μ))) )dΩᶜ

  # Rigid body motion and volum constraint
  RB¹ = VectorValue(1.0,0.0)
  r¹(u,ℓ) = ∫( ( RB¹⋅(ℓ*u) )*y )dΓ
  r²(u,ℓ) = ∫( ( u⋅(ℓ*nΓ ) )*y )dΓ

  aᵛ((υ,l¹,l²),(μ,ℓ¹,ℓ²)) = 
    aʷ(υ,μ,eₕ) + aᶠ(υ,μ,R) + sᵘ(υ,μ) + r¹(υ,ℓ¹) + r¹(μ,l¹) + r²(υ,ℓ²) + r²(μ,l²)
  bᵛ((μ,ℓ¹,ℓ²)) = f(μ, ρₕ,R,eₕ)

  aᵛ, bᵛ 
end 

function turnover_axisymmetric(u,eₕ,dΓ,dΩᶜ,nΓ,De,
    dt::Float64,γ::Float64,τᵈkₒ::Float64)
  
  m(e,ε)  = ∫( (1/dt)*(e*ε)*y )dΓ
  sᵈ(e,ε) = ∫( ( ∇ᵈ(e,nΓ)⋅∇ᵈ(ε,nΓ) )*y )dΓ
  c(e,ε)  = ∫( ( (u⋅∇ᵈ(e,nΓ))*ε + (tr(∇ᵈ(u,nΓ))+u⋅iy)*(e*ε) )*y )dΓ #TODO Ask Eric about this
  #advection(rac,w,v) =  ∫( (w*( ∇ᵈ(rac,nΓ)⋅v + rac*divᶜ(v,nΓ))  )*y )dΓ
  r(e,ε)  = ∫( τᵈkₒ*(e*ε)*y )dΓ
  d(e,ε) = ∫( ( De * (∇ᵈ(e,nΓ)⋅∇ᵈ(ε,nΓ)) )*y )dΓ
  l(ε)    = ∫( τᵈkₒ*ε*y )dΓ
#  s(υ,μ)  = ∫( γ*((nΓ⋅∇(υ))⊙(nΓ⋅∇(μ))) )dΩᶜ

  aᵉ(e,ε) = m(e,ε) + c(e,ε) + r(e,ε) + d(e,ε) + sᵈ(e,ε) #+ s(e,ε)
  bᵉ(ε)   = m(eₕ,ε) + l(ε)

  aᵉ,bᵉ
end

function cortical_flow_problem_axisymmetric(
    ρₕ,R,dΩᶜ,dΓ,nΓ,γ::Float64,Pe::Float64,χᵣ::Float64,χ₀::Float64,ξ₀,σₐ⁰,sigmaₐ⁰,  sigmaρ⁰, sigmaR⁰)

  # Viscous term
  aʷ(u,v) = 
    ∫( ( εᶜ(u,nΓ)⊙εᵈ(v,nΓ) + divᶜ(u,nΓ)⋅divᶜ(v,nΓ) + 
         2*(u⋅iy)*(v⋅iy) + divᶜ(u,nΓ)*(v⋅iy) + 
         divᶜ(v,nΓ)*(u⋅iy) )*y )dΓ

  χ(R) = (χ₀+χᵣ*R)

  # Friction term
  aᶠ(u,v,R) = ∫(χ(R)*(u⋅v)*y )dΓ

  # Activity function (relates Rho concentration to active stress)
  #ξ(e) = 2.0 * e*e / ( 1.0 + e*e )
  #ξ(ρ) = ∇ᵈ(ρ,nΓ)#⋅(TensorValue(0.0,-1.0,1.0,0.0)⋅nΓ)#divᶜ( ρ,nΓ)
  #TODO  
  function sigmaₐ(ρ,R)
      sigmaₐ = sigmaₐ⁰ + sigmaρ⁰ * ρ - sigmaR⁰ * R
      sigmaₐ > 0 ? sigmaₐ : zero(typeof(sigmaₐ))
  end
  # Active force term
  f(μ, ρ,R) = ∫( ( -(divᶜ(μ,nΓ)+μ⋅iy)*(sigmaₐ∘(ρ,R))*σₐ⁰ ) * ξ₀ )dΓ
# f(μ,e) = ∫( Pe * ( -(divᶜ(μ,nΓ)+μ⋅iy)*(ξ∘(e)) ) * ξ₀ )dΓ OG FUNCTION
# f(μ,e) = ∫( Pe * ( -(divᶜ(μ,nΓ))*(ξ∘(e)) ) * ξ₀ )dΓ TRANSPORT VERSION

  # Stabilisation term for velocity
  sᵘ(υ,μ) = ∫( γ * ((nΓ⋅ε(υ))⊙(nΓ⋅ε(μ))) )dΩᶜ

  # Rigid body motion and volum constraint
  RB¹ = VectorValue(1.0,0.0)
  r¹(u,ℓ) = ∫( ( RB¹⋅(ℓ*u) )*y )dΓ
  r²(u,ℓ) = ∫( ( u⋅(ℓ*nΓ ) )*y )dΓ

  aᵛ((υ,l¹,l²),(μ,ℓ¹,ℓ²)) =
    aʷ(υ,μ) + aᶠ(υ,μ,R) + sᵘ(υ,μ) + 
    r¹(υ,ℓ¹) + r¹(μ,l¹) + r²(υ,ℓ²) + r²(μ,l²)
  bᵛ((μ,ℓ¹,ℓ²)) = f(μ, ρₕ,R)

  aᵛ, bᵛ
end


function cortical_flow_problem_mechanochemical_axisymmetric(
    ρₕ,ez,dΩᶜ,dΓ,nΓ,γ::Float64,Pe::Float64,χᵣ::Float64,χ₀::Float64,ξ₀,sigmaₐ⁰,  sigmaρ⁰)

  # Viscous term
  aʷ(u,v) = 
    ∫( ( εᶜ(u,nΓ)⊙εᵈ(v,nΓ) + divᶜ(u,nΓ)⋅divᶜ(v,nΓ) + 
         2*(u⋅iy)*(v⋅iy) + divᶜ(u,nΓ)*(v⋅iy) + 
         divᶜ(v,nΓ)*(u⋅iy) )*y )dΓ

  χ(R) = (χ₀+χᵣ*R)

  # Friction term
  aᶠ(u,v,R) = ∫(χ(R)*(u⋅v)*y )dΓ

  # Activity function (relates Rho concentration to active stress)  
  function sigmaₐ(ρ,R)
      sigmaₐ = sigmaₐ⁰ + sigmaρ⁰ * ρ #- sigmaR⁰ * R
      sigmaₐ > 0 ? sigmaₐ : zero(typeof(sigmaₐ))
  end
  # Active force term
  f(μ, ρ,R) = ∫( ( -(divᶜ(μ,nΓ)+μ⋅iy)*(sigmaₐ∘(ρ,R)) ) * ξ₀ )dΓ 

  # Stabilisation term for velocity
  sᵘ(υ,μ) = ∫( γ * ((nΓ⋅ε(υ))⊙(nΓ⋅ε(μ))) )dΩᶜ

  # ** weak tangentiality **
  η = 10.0 / ((2/40)^2)
  k(u,v) = ∫( η*((u⋅nΓ)*(v⋅nΓ)) )dΓ

  # Rigid body motion and volum constraint
  RB¹ = VectorValue(1.0,0.0)
  r¹(u,ℓ) = ∫( ( RB¹⋅(ℓ*u) )*y )dΓ
  r²(u,ℓ) = ∫( ( u⋅(ℓ*nΓ ) )*y )dΓ

  aᵛ(υ,μ) = aʷ(υ,μ) + sᵘ(υ,μ) + k(υ,μ) + aᶠ(υ,μ,ez)
  bᵛ(μ) = f(μ, ρₕ,ez)

  aᵛ, bᵛ
end

function cortical_flow_problem_axisymmetric(ulₕ,plₕ,eₕ,dΩᶜ,dΓ,nΓ,
    γ::Float64,Pe::Float64,μˡ::Float64,R::Float64,ξ₀::Function)

  aʷ(υ,μ) = ∫( 2.0 * ( εᶜ(υ,nΓ)⊙εᵈ(μ,nΓ) + (υ⋅iy)*(μ⋅iy) ) * y )dΓ

  ξ(e) = 2.0 * e*e / ( 1.0 + e*e )

  f(μ,e) = ∫( Pe * ( -(divᶜ(μ,nΓ)+μ⋅iy)*(ξ∘(e)) ) * ξ₀ )dΓ

  σᵘ(ε,q) = 2.0 * μˡ * R * ε - q * one(ε)
  βʳ(μ,u,p) = ∫( ( μ⋅((σᵘ∘(ε(u),p))⋅nΓ) ) * y )dΓ

  sᵘ(υ,μ) = ∫( γ * ((nΓ⋅ε(υ))⊙(nΓ⋅ε(μ))) )dΩᶜ

  RB¹ = VectorValue(1.0,0.0)

  r¹(u,ℓ) = ∫( ( u⋅(ℓ*RB¹) )*y )dΓ
  r²(u,ℓ) = ∫( ( u⋅(ℓ*nΓ ) )*y )dΓ

  aᵛ((υ,l¹,l²),(μ,ℓ¹,ℓ²)) =
    aʷ(υ,μ) + sᵘ(υ,μ) + r¹(υ,ℓ¹) + r¹(μ,l¹) + r²(υ,ℓ²) + r²(μ,l²)
  bᵛ((μ,ℓ¹,ℓ²)) = f(μ,eₕ) - βʳ(μ,ulₕ,plₕ)

  aᵛ, bᵛ
end

function bulk_flow_problem_axisymmetric(
    υ,dΩ,dΓ,nΓ,μˡ::Float64,R::Float64,γ::Float64,h::Float64)
  
  aᵇ(u,v) = ∫( 2.0 * μˡ * R * (ε(u)⊙ε(v) + (u⋅iy)*(v⋅iy)) * y )dΩ

  bᵇ(v,q) = ∫( q*(∇⋅v + v⋅iy)*y )dΩ

  σᵘ(ε,q) = 2.0 * μˡ * R * ε - q * one(ε)
  αˡ(u,v,p,q) = ∫( ( (γ/h)*(u⋅v)           -
                      u⋅((σᵘ∘(ε(v),q))⋅nΓ) -
                      v⋅((σᵘ∘(ε(u),p))⋅nΓ) )*y )dΓ
  αʳ(v,q,υ)   = ∫( ( (γ/h)*(υ⋅v) - υ⋅((σᵘ∘(ε(v),q))⋅nΓ) )*y )dΓ

  r(p,ℓ) = ∫( ( p*ℓ ) * y )dΩ

  aᵘ((uˡ,pˡ,l),(vˡ,qˡ,ℓ)) = 
    aᵇ(uˡ,vˡ) - bᵇ(vˡ,pˡ) - bᵇ(uˡ,qˡ) + 
    αˡ(uˡ,vˡ,pˡ,qˡ) + r(pˡ,ℓ) + r(qˡ,l)
  bᵘ((vˡ,qˡ,l)) = αʳ(vˡ,qˡ,υ)

  aᵘ, bᵘ
end

function transport_problem_axisymmetric(u,eₕ,dΓ,dΩᶜ,nΓ,
    dt::Float64,γ::Float64,τᵈkₒ::Float64)
  
  m(e,ε)  = ∫( (1/dt)*(e*ε)*y )dΓ
  sᵈ(e,ε) = ∫( ( ∇ᵈ(e,nΓ)⋅∇ᵈ(ε,nΓ) )*y )dΓ
  c(e,ε)  = ∫( ( (u⋅∇ᵈ(e,nΓ))*ε + (tr(∇ᵈ(u,nΓ))+u⋅iy)*(e*ε) )*y )dΓ
  r(e,ε)  = ∫( τᵈkₒ*(e*ε)*y )dΓ
  l(ε)    = ∫( τᵈkₒ*ε*y )dΓ
  s(υ,μ)  = ∫( γ*((nΓ⋅∇(υ))⊙(nΓ⋅∇(μ))) )dΩᶜ

  aᵉ(e,ε) = m(e,ε) + c(e,ε) + r(e,ε) + sᵈ(e,ε) + s(e,ε)
  bᵉ(ε)   = m(eₕ,ε) + l(ε)

  aᵉ,bᵉ
end

function cortical_flow_problem_3D(ulₕ,plₕ,eₕ,dΩᶜ,dΓ,nΓ,
    γ::Float64,Pe::Float64,μˡ::Float64,R::Float64,ξ₀::Function)

  β = 0.001
  aʷ(υ,μ) = ∫( 2.0 * ( εᶜ(υ,nΓ)⊙εᵈ(μ,nΓ) ) + β * (υ⋅μ) )dΓ

  ξ(e) = 2.0 * e*e / ( 1.0 + e*e )

  f(μ,e) = ∫( Pe * ( -(divᶜ(μ,nΓ))*(ξ∘(e)) ) * ξ₀ )dΓ

  σᵘ(ε,q) = 2.0 * μˡ * R * ε - q * one(ε)
  βʳ(μ,u,p) = ∫( μ⋅((σᵘ∘(ε(u),p))⋅nΓ) )dΓ

  sᵘ(υ,μ) = ∫( γ*((nΓ⋅ε(υ))⊙(nΓ⋅ε(μ))) )dΩᶜ

  RB¹ = VectorValue(1.0,0.0,0.0)
  RB² = VectorValue(0.0,1.0,0.0)
  RB³ = VectorValue(0.0,0.0,1.0)
  RB⁴ = x -> VectorValue(0.0,-x[3],x[2])
  RB⁵ = x -> VectorValue(x[3],0.0,-x[1])
  RB⁶ = x -> VectorValue(-x[2],x[1],0.0)

  r¹(u,ℓ) = ∫( (u⋅RB¹)*ℓ )dΓ
  r²(u,ℓ) = ∫( (u⋅RB²)*ℓ )dΓ
  r³(u,ℓ) = ∫( (u⋅RB³)*ℓ )dΓ
  r⁴(u,ℓ) = ∫( (u⋅RB⁴)*ℓ )dΓ
  r⁵(u,ℓ) = ∫( (u⋅RB⁵)*ℓ )dΓ
  r⁶(u,ℓ) = ∫( (u⋅RB⁶)*ℓ )dΓ
  r⁷(u,ℓ) = ∫( (u⋅nΓ )*ℓ )dΓ

  aᵛ((υ,l¹,l²,l³,l⁴,l⁵,l⁶,l⁷),(μ,ℓ¹,ℓ²,ℓ³,ℓ⁴,ℓ⁵,ℓ⁶,ℓ⁷)) =
    aʷ(υ,μ) + sᵘ(υ,μ) + 
    r¹(υ,ℓ¹) + r¹(μ,l¹) + r²(υ,ℓ²) + r²(μ,l²) +
    r³(υ,ℓ³) + r³(μ,l³) + r⁴(υ,ℓ⁴) + r⁴(μ,l⁴) +
    r⁵(υ,ℓ⁵) + r⁵(μ,l⁵) + r⁶(υ,ℓ⁶) + r⁶(μ,l⁶) +
    r⁷(υ,ℓ⁷) + r⁷(μ,l⁷)
  bᵛ((μ,ℓ¹,ℓ²,ℓ³,ℓ⁴,ℓ⁵,ℓ⁶,ℓ⁷)) = f(μ,eₕ) - βʳ(μ,ulₕ,plₕ)

  aᵛ, bᵛ
end

function bulk_flow_problem_3D(
  υ,dΩ,dΓ,nΓ,μˡ::Float64,R::Float64,γ::Float64,h::Float64)

  aᵇ(u,v) = ∫( 2.0 * μˡ * R * (ε(u)⊙ε(v)) )dΩ

  bᵇ(v,q) = ∫( q*(∇⋅v) )dΩ

  σᵘ(ε,q) = 2.0 * μˡ * R * ε - q * one(ε)
  αˡ(u,v,p,q) = ∫( (γ/h)*(u⋅v)           -
                    u⋅((σᵘ∘(ε(v),q))⋅nΓ) -
                    v⋅((σᵘ∘(ε(u),p))⋅nΓ) )dΓ
  αʳ(v,q,υ) = ∫( (γ/h)*(υ⋅v) - υ⋅((σᵘ∘(ε(v),q))⋅nΓ) )dΓ

  r(p,ℓ) = ∫( p*ℓ )dΩ

  aᵘ((uˡ,pˡ,l),(vˡ,qˡ,ℓ)) = 
    aᵇ(uˡ,vˡ) - bᵇ(vˡ,pˡ) - bᵇ(uˡ,qˡ) + 
    αˡ(uˡ,vˡ,pˡ,qˡ) + r(pˡ,ℓ) + r(qˡ,l)
  bᵘ((vˡ,qˡ,l)) = αʳ(vˡ,qˡ,υ)

  aᵘ, bᵘ
end

function transport_problem_3D(u,eₕ,dΓ,dΩᶜ,nΓ,
    dt::Float64,γ::Float64,τᵈkₒ::Float64)

  m(e,ε)  = ∫( (1/dt)*(e*ε) )dΓ
  sᵈ(e,ε) = ∫( ∇ᵈ(e,nΓ)⋅∇ᵈ(ε,nΓ) )dΓ
  c(e,ε)  = ∫( (u⋅∇ᵈ(e,nΓ))*ε + tr(∇ᵈ(u,nΓ))*(e*ε) )dΓ
  r(e,ε)  = ∫( τᵈkₒ*(e*ε) )dΓ
  l(ε)    = ∫( τᵈkₒ*ε )dΓ
  s(υ,μ)  = ∫( γ*((nΓ⋅∇(υ))⊙(nΓ⋅∇(μ))) )dΩᶜ

  aᵉ(e,ε) = m(e,ε) + c(e,ε) + r(e,ε) + sᵈ(e,ε) + s(e,ε)
  bᵉ(ε)   = m(eₕ,ε) + l(ε)

  aᵉ,bᵉ
end