# Standard code to run a simulation using the mechanochemical model presented in article "INSERT DOI"
# The system will start with flat Rho and Rac at an steady state
# Different mechanical parametes can eb changed in the main section of the code.
# One can work with pure local inhibition by changing the following mechanical parameters:
# vCTE=0   α=0   β=0   σₐ₀=0
# Written by Andreu F Gallen working in Turlier lab and in collaboration with Orion Weiner's lab

using Gridap
using GridapEmbedded
using Plots
using Plots.PlotMeasures
using SurfaceBulkViscousFlows


###############CONSTANT MECHANICAL PARAMETERS ##############
partition=40       #size of the FE system
χ = 100              # [pN s/um^3] friction coefficient
η=1000              # [pN s/ um] 2D viscosity of the cortex
L = 12 #1.2
L2 = 12 #1.2
R = 10       # [um] System size length
R2 = 10      # [um] System size length
k=100.0            # [pN /um] elastic constant for the membrane
σₐ₀ =  20.0        # [pN /um^2] maximum active "pressure"

koff = 1.4  # ezrin koff [1/s] toff=0.7s Fritzsche et al
kon  = 5.0  # ezrin kon [1/s] ton =0.2s Fritzsche et al
D = 0.3   # diffusion of the membrane [um^2/s] 0.003  Fritzsche et al
M0 = 0.1 #total amount of ezrin
const wrac=R2*π/2.5 #Opto signal use this for the Gaussian distribution width
#in case one wants nonlinear deactication if the system would be unstable
λᵇ = 0.001   # for bound ezrin 
λʳᴬ= 0   #for Rho

#adimensional numbers that define the system
τₐ = η/σₐ₀ # 100s for σₐ₀=100
τₑ = η/k #almost always?
Pe⁻¹ = D*η/(σₐ₀*R^2) #D*η/(σₐ₀*L^2) 0.00075 for σₐ₀=100
λ = sqrt( η*R*π /(χ*R*R*π*π*M0) ) #\bar λ in the suppl material. Adimensional lenthscale

tenth=2.0 #threshold for tension to activate rho
sig0 = 2.00  #adimensionalizes the tension term for rho
MCAbth = 0.01 #MCA/ezrin value below which protrusion starts
rho0 = 0.01 #adimensionalizes the MCA concentration term for rac

#Rac Coefficients
dᵃ = 0.04  #deactivation
α₀ = 3.00  #basal activation
α =0.2*α₀ #activation through mechanics, zero for only local inhibition
Drac = 0.005 #diffusion rate
αopto = 20 # opto input for Rac
a_t = 0

#Rho Coefficients
dᵇ = 0.04 #deactivation
β₀ = 3.00 #basal activation
β = 0.2*β₀ #activation through mechanics, zero for only local inhibition
Drho = 0.005 #diffusion rate 
βopto = 0 # opto input for Rho
b_t = 0

Dₑ = 0.01
kₑoff = 1.0
kₑon = 1.0 

#time variables
topto=20 #time to start opto signal
Δt  = 1.0 #timestep
T = 400 #time to finish simulation
#if we wanna do multiple simulations in a row with different variables
setsv=1
setsx=2 

#output folder name
len=trunc(λ,digits=2)
 @info "Characteristic length $len system size $R"
simulation = "singlet/T=$T part=$partition db=$dᵇ da=$dᵃ a0=$α₀ b0=$β₀ bopto=$βopto aopto=$αopto a_t=b_t=$b_t/sig_a=$σₐ₀ len=$len rho0=$rho0 D=$D Drac=Drho=$Drac M0=$M0 deltat=$Δt/"
#store VTUs in one folder
pVTU="./VTU/"*simulation
mkpath(pVTU)
#Store PNGs in another
pPNG="./PNG/"*simulation
mkpath(pPNG)

#IF we run multiple simulations with changing parameters, we use this arrays to store the results of each simulation
xχ = zeros(setsx)
xη = zeros(setsv)
tension = zeros(setsx,setsv,Int32(T/Δt)+1,partition-2)
v = zeros(setsx,setsv,Int32(T/Δt)+1,partition-1)
MCA_b = zeros(setsx,setsv,Int32(T/Δt)+1,partition+1)
a = zeros(setsx,setsv,Int32(T/Δt)+1,partition+1)
b = zeros(setsx,setsv,Int32(T/Δt)+1,partition+1)
tension2 = zeros(setsx,setsv,Int32(T/Δt)+1,partition-2)
v2 = zeros(setsx,setsv,Int32(T/Δt)+1,partition-1)
MCAb2 = zeros(setsx,setsv,Int32(T/Δt)+1,partition+1)
αt = zeros(setsx,setsv,Int32(T/Δt)+1,partition)
βt = zeros(setsx,setsv,Int32(T/Δt)+1,partition)

#FOR for the possible multiple simulations
for i in 2:1:setsx
  for j in 1:1:setsv
    # χ = 2*10^((i))    # [pN s/um^3] friction coefficient
    # χ = trunc(χ,digits=2)
    λ⁻² = 1/(λ*λ)
    xχ[i]=χ 
    #η= 0+0*j+10^(1*(j+1))              # [pN s/ um] 2D viscosity of the cortex
    lη = η  
    xη[j]=lη
  #  MCA_b[i,j,:,:],v[i,j,:,:],a[i,j,:,:],b[i,j,:,:]= 
run_singlet_axisymmetric(χ,lη,T,Δt,partition,
        L,simulation,wrac,αopto,βopto,kon,koff,M0,α₀,β₀,
        σₐ₀,λᵇ,D,Drac,Drho,a_t,b_t,α,β,dᵃ,dᵇ,
        sig0,tenth,MCAbth,topto,R,R2,L2)
  end
end

