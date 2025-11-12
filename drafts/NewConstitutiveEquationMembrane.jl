# HOMEWORK
# 
# Andreu: 1. Checks calculations by Eric on notes
#         2. Finishes implementation of the weak form 
#            of the new membrane equation.
# Eric: 1. Checks relation between surface gradient and ∂x/∂θ
#       2. Systematic linearisation with directional derivative

# Weak form of the membrane equation
#
# Nonlinear strain rate: 
# εᴾ(u) = ε(u) + 0.5 * ∇(u)ᵗ⋅∇(u) = εᴾ(u) = ε(u) + εᴺ(u)
#
# Terms of the bilinear and linear forms
# TERM 1. ∫( 2M⋅εᴾ(u):ε(v) )dΓ = ∫( 2M⋅ε(u):ε(v) + 
#                                   2M⋅εᴺ(u):ε(v) )dΓ
#
# OBS 1. Eric checks relation ∇ᵈ(x,nΓ) and ∂x/∂θ
# OBS 2. To implement function θ
# OBS 3. Beware of orientation of θ, assuming 
#        that θ = 0 at the North Pole, 
#        and θ = π at the South Pole.
#        > Implement θ such that θ = 0 at right Pole
#          and θ = π at left Pole.
# OBS 4. Terms like x*x*w can be linearised as
#        x_old*x*w or x_old^2*w. Eric will check
#        how to rigorously linearise these terms.
aᴹ(M,R,x_old,x,w,θ) = 
  ∫( ( 2*M * ( x*w/2 + 
               R^2 * ( ∇ᵈ(x,nΓ)⋅∇ᵈ(w,nΓ) ) + 
               (cot∘(θ)^2) * x*w ) ) * sin∘(θ) )dΓ +
  ∫( ( M/R * ( 
    ( x_old * x + R^2 * ( ∇ᵈ(x_old,nΓ)⋅∇ᵈ(x,nΓ) ) ) * ( R * ∇ᵈ(w,nΓ) ) + 
    ( cot∘(θ)^3 * x_old ) * x * w ) ) * sin∘(θ) )dΓ
#
# TERM 2. ∫( L⋅(tr(εᴾ(u))Id):ε(v) )dΓ = ∫( L⋅tr(ε(u)):ε(v) +  
#                                          L⋅tr(εᴺ(u)):ε(v) )dΓ
# Homework: Implement TERM 2
aᴸ(L,R,x_old,x,w,θ) = 
  ∫( L*(R^2*(∇ᵈ(x,nΓ)⋅∇ᵈ(w,nΓ)) + cot∘(θ)*x*∇ᵈ(w,nΓ) + cot∘(θ)*v*∇ᵈ(x,nΓ) + 
               (cot∘(θ)^2) * x*w  ) * sin∘(θ) )dΓ +
  ∫( 0.5/R*L*( ( (1+(cot∘(θ))^2)*x_old*x + R*∇ᵈ(x,nΓ)*R*∇ᵈ(x_old,nΓ) )*( R*∇ᵈ(w,nΓ) 
      + cot∘(θ)*v ) ) * sin∘(θ) )dΓ
#cosecant(θ) = 1.0 / sin(θ)
#
# TERM 3. ∫( div⋅(S⋅ID) )dΓ = 0 because S real constant
# However, add basal stress S⋅ID when reporting stresses
#

# New membrane equation
a(L,M,R,x_old,x,w,θ) = aᴸ(L,R,x_old,x,w,θ) + aᴹ(M,R,x_old,x,w,θ)

# Preserve mass term for Backward Euler time integration
m(MCA_b,Δt,x,w) = ∫( ( (χ*MCA_b) * (x*w) / Δt )*y )dΓ

A(x,w) = m(uh_MCAb,Δt,x,w) + a(L,M,R,x_old,x,w,θ) + s₀x(x,w)
B(w) = m(uh_MCAb,Δt,uh_x,w) + bₓ(uh_MCAb,uh_v,w) 
op_x = AffineFEOperator(A,B,X,WD0)
uh_x = solve(op_x)