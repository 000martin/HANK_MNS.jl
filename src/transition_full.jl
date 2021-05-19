#features files used to compute the transition path of the full model

using Plots

"""
    get_transition_full()

Solves for a steady state and a transition path.
Inputs: TR (time until single-quarter interest rate change), T (time horizon transition path), RChange (Change of R at time TR)
"""   
function get_transition_full(TR::Int,T::Int,p::params,SS::steady_state; RChange::Float64 = -0.005)

    @unpack Rbar = p

    #generate R path
    Rpath = repeat([Rbar],T); Rpath[TR+1] = Rbar + RChange
    
    #make guess for paths: just steady state state_values for entire transition periods
    wpath = repeat([SS.w],T) ; div_path = repeat([SS.div],T) 
    Spath = ones(T) ; 

    #get transition path
    tp = solve_for_transition(Rpath,wpath,div_path,Spath,SS,p)

    return tp 
end

"""
    solve_for_transition(Rpath::Array{Float64,1},wpath::Array{Float64,1},div_path::Array{Float64,1},Spath::Array{Float64,1},SS::steady_state,p::params)

Solves for the perfect foresight transition path for a given interest rate path. wpath, divpath, etc. are guesses for the transition path
Note that guesses should include SS values as final element.
"""

function solve_for_transition(Rpath::Array{Float64,1},wguess::Array{Float64,1},div_path::Array{Float64,1},
                            Spath::Array{Float64,1},SS::steady_state,p::params; S_tol::Float64= 1e-6, w_tol::Float64 = 1e-6)

 #unpack some parameters
 @unpack B,tax_weights, Γ, ψ, μ, β , θ , nk, nz, nb = p


 #back out number of periods
 T = length(Rpath)   

 #get path for tax rate
 τ_path = (B*SS.Y/(Γ'*tax_weights))*(1.0 .- 1.0./Rpath)
 wpath = wguess[:] ; 

 #make objects exist outside while loop
 Ypath = ones(T-1) ; pΠpath = ones(T)  
 Dpath = Array{Float64,2}(undef,nk*nz,T) ; cpol_path =  Array{Float64,2}(undef,nb*nz,T)

 #outer loop solving for S (price dispersion) path
 iterS = 1 ; distS = 1.0
 while (iterS < 100) & (distS > S_tol )

    
    #inner loop solving for wage path
    iterW = 1 ; distW = 1.0
    while (iterW < 100) & (distW > w_tol) #does not need to converge fully for every S iteration
    
    #solve HH problem backwards
    cs = solveback(reshape_c(SS.c_policies,p),wpath,Rpath,τ_path,div_path,β,p)
    cpol_path = cs

    #use result to simulate aggregate forwards
    Cpath,Lpath,Bpath = simulate_forward(SS.D,cpol_path,Rpath,wpath,div_path,τ_path,p)

    #define distribution path
    #Dpath = Ds

    #define output path
    Ypath = Cpath[:]

    #get aggregate labor demand path
    Npath = Spath[1:T-1].*Ypath

    #update wage and dividend paths
    oldwage = wpath[:]
    wpath[1:T-1] = wpath[1:T-1].*(Npath[1:T-1]./Lpath[1:T-1]).^ψ
    
    #same update rule as in MNS
    wpath[1:T-1] = 0.25*wpath[1:T-1] .+ 0.75*oldwage[1:T-1]

    div_path[1:T-1] = Ypath[1:T-1] .- wpath[1:T-1].*Npath[1:T-1]

    #calculate distance
    distW = maximum(abs.(wpath[1:T-1]./oldwage[1:T-1] .- 1.0))

    println("Current W distance: ", distW," Current wage iteration: ",iterW)
    iterW = iterW + 1
    end

    #initialize pbarA and pbarB terms (nominator and denominator in MNS eqaution (7))
    pbarA = μ*SS.w*SS.Y / (1-β*(1-θ)) ; pbarB = SS.Y/(1-β*(1-θ))

    #pre-allocate
    pΠpath = ones(T) ; pstar = ones(T-1) 

    #solve backwards for pbarA, pbarB and reset inflation
    for t = T-1:-1:1
        pbarA = μ*wpath[t]*Ypath[t] + β*(1-θ)*(pΠpath[t+1]^(μ/(μ-1)))*pbarA
        pbarB = Ypath[t] + β*(1-θ)*(pΠpath[t+1]^(1/(μ-1)))*pbarB
        pstar[t] = pbarA/pbarB
        pΠpath[t] = ((1-θ)/(1-θ*pstar[t]^(1/(1-μ))))^(1-μ)
    end

    oldS = Spath ; Spath = ones(T)
    Slast = 1.0 #steady state price dispersion
    #solve for S path
    for t = 1:T-1
        Spath[t] = (1-θ)*Slast*pΠpath[t]^(μ/(μ-1)) + θ*pstar[t]^(μ/(1-μ))
        Slast = Spath[t]
    end

    #compute Distance
    distS = maximum(abs.(Spath./oldS .- 1.0))
    println("Current S distance: ", distS," Current S iteration: ",iterS)
    println(" ")

    iterS = iterS + 1
 end

 return transition_full(Spath,wpath,pΠpath,Ypath,Rpath,τ_path,div_path)

end



"""
    simulate_forward()

Simulates a distribution of households and computes total consumption (C), labor supply (L), asset position (B) and wealth distribution (D)
"""

function simulate_forward(D0::Array{Float64,1},cpol_path::Array{Float64,2},Rpath::Array{Float64,1},
                        wpath::Array{Float64,1},div_path::Array{Float64,1},τ_path::Array{Float64,1},p::params)

 @unpack k_grid = p

 #back out number of periods
 T = length(Rpath)

 #pre-allocate some arrays
 Cpath = Array{Float64,1}(undef,T-1);  Lpath = Array{Float64,1}(undef,T-1)
 Bpath = Array{Float64,1}(undef,T-1)  
 Dpath = Array{Float64,2}(undef,p.nz*p.nk,T)
 Dpath[:,1] .= D0 

 for t = 1:T-1

 #simulate forward
  cpols = reshape_c(cpol_path[:,t],p)
  Cpath[t],Lpath[t] = aggregate_C_L(Dpath[:,t],cpols,Rpath[t],wpath[t],τ_path[t],div_path[t],p)
  Dpath[:,t+1]     .= forward_dist(Dpath[:,t],forwardmat(cpols,Rpath[t],wpath[t],τ_path[t],div_path[t],p))
  Bpath[t]          = dot(Dpath[:,t+1],repeat(k_grid,3))

 #check distribution
 @assert (abs(sum(Dpath[:,t+1]) - 1.0) < 1e-6)
 end
 

 return Cpath , Lpath, Bpath

    
end

"""
    simulate_step()

Conducts forward simulation for one period. Helper function to simulate_forward, equivalent to simulatestep() in MNS code.
Originally used in loop in Simulate_forward, replaced it to do pre-allocation there.

The following was originally in the simulate forward loop
# Cpath[t], Lpath[t], Bpath[t], Dpath[:,t+1]  = simulate_CLB(Dpath[:,t],reshape_c(cpol_path[:,t],p)
#                                                  ,Rpath[t],wpath[t],τ_path[t],div_path[t],p)

"""
function simulate_step(D::Array{Float64,1},c_pol::Array{Float64,2},R::Float64,w::Float64,τ::Float64,div::Float64,p::params)

    @unpack k_grid = p

 #household consumption and labor supply
 C,L = aggregate_C_L(D,c_pol,R,w,τ,div,p)

 Pi = forwardmat(c_pol,R,w,τ,div,p)

 #distribution in next period
 Dprime = Pi'*D

 #aggregate assets
 Assets = dot(Dprime,repeat(k_grid,3))

 #test for validity of distribution 
 @assert (abs(sum(Dprime) - 1.0) < 1e-6)

 return C, L, Assets, Dprime
    
end

"""
    forward_dist

Calculates Asset distribution in next period given transition matrix Pi and current distribution D
"""
function forward_dist(D::Array{Float64,1},Pi::SparseMatrixCSC)
return Pi'*D
end