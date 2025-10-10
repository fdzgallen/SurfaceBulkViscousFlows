iy(x) = VectorValue( 0.0, 1.0 / x[2] ); y(x) = x[2]


function rac_rho_weak_forms2(Δt,dᵃ,dᵇ,Drac,Drho,nΓ,dΓ)
  #DEFINING the equations for Rac and Rho
  m2(Δt,A,w) = ∫( ( (A*w)/Δt )*y )dΓ
  a_rac(rac,w) = (1/dᵃ) * m2(Δt,rac,w) + ∫( ( Drac * (∇ᵈ(rac,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ + ∫( ( w*rac )*y )dΓ
  b_rac(w,rho,α₀v,rac_old) = (1/dᵃ) * m2(Δt,rac_old,w) +   ∫( ( w*(α₀v)/(1+rho*rho) )*y )dΓ   
 
  a_rho(rho,w) = (1/dᵇ)*m2(Δt,rho,w) + 
    ∫( ( Drho * (∇ᵈ(rho,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ + ∫( ( w*rho )*y )dΓ # + ∫( ( λʳᴬ*((rho*rho*rho)*w) )*y )dΓ  
  b_rho(w,rac,β₀v,rho_old) = (1/dᵇ)*m2(Δt,rho_old,w) + ∫( ( w*(β₀v)/(1+rac*rac) )*y )dΓ 

  a_rac, b_rac, a_rho, b_rho
end


function rac_rho_weak_forms_conserved(Δt,dᵃ,dᵇ,Drac,Drho,nΓ,dΓ)
  #DEFINING the equations for Rac and Rho
  #TODO add y's, convert into axisymmetric
  m2(Δt,A,w) = ∫( ( (A*w)/Δt )*y )dΓ
  a_rac(rac,w) = (1/dᵃ)*m2(Δt,rac,w) + ∫( ( Drac * (∇ᵈ(rac,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ + ∫(w*rac)dΓ  
  b_rac(w,rho,α₀v,rac_i,rac_old) =  (1/dᵃ)*m2(Δt,rac_old,w) + ∫( w*(rac_i)*(α₀v/(1+rho*rho))*y )dΓ  
  a_rac_i(rac_i,w) = ∫(((w*(rac_i))*y))dΓ  
  b_rac_i(w,a_sum) = ∫((w*(a_t - a_sum))*y)dΓ

  a_rho(rho,w) = (1/dᵇ)*m2(Δt,rho,w) + ∫( ( Drho * (∇ᵈ(rho,nΓ)⋅∇ᵈ(w,nΓ)) )*y )dΓ  + 
  ∫(w*rho)dΓ  #  +∫(λʳᴬ*((rho*rho*rho)*w))dΓ  
  b_rho(w,rac,rho_i) =  (1/dᵇ)*m2(Δt,uh_rho_old,w) + ∫( (w*(rho_i)*(β₀v/(1+rac*rac)))*y)dΓ  
  a_rho_i(rho_i,w) = ∫((w*(rho_i))*y)dΓ 
  b_rho_i(w,b_sum) = ∫((w*(b_t - b_sum))*y)dΓ

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

function cortical_flow_problem_axisymmetric(
    ρₕ,R,dΩᶜ,dΓ,nΓ,γ::Float64,Pe::Float64,χ::Float64,χ₀::Float64,ξ₀)

  # Viscous term
  aʷ(u,v) = 
    ∫( ( εᶜ(u,nΓ)⊙εᵈ(v,nΓ) + divᶜ(u,nΓ)⋅divᶜ(v,nΓ) + 
         2*(u⋅iy)*(v⋅iy) + divᶜ(u,nΓ)*(v⋅iy) + 
         divᶜ(v,nΓ)*(u⋅iy) )*y )dΓ

  # Friction term
  aᶠ(u,v,R) = ∫( (χ₀+χ*R)*(u⋅v)*y )dΓ

  # Activity function (relates Rho concentration to active stress)
  #ξ(e) = 2.0 * e*e / ( 1.0 + e*e )
  #ξ(ρ) = ∇ᵈ(ρ,nΓ)#⋅(TensorValue(0.0,-1.0,1.0,0.0)⋅nΓ)#divᶜ( ρ,nΓ)

  # Active force term
  f(μ, ρ,R) = ∫( ( -(divᶜ(μ,nΓ)+μ⋅iy)*(∇ᵈ(ρ,nΓ)⋅(TensorValue(0.0,-1.0,1.0,0.0)⋅nΓ) - ∇ᵈ(R,nΓ)⋅(TensorValue(0.0,-1.0,1.0,0.0)⋅nΓ)) ) * ξ₀ )dΓ

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