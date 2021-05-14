#features files used to compute the transition path of the full model

"""
    solve_for_transition(Rpath::Array{Float64,1},wpath::Array{Float64,1},div_path::Array{Float64,1},Spath::Array{Float64,1},SS::steady_state,p::params)

Solves for the perfect foresight transition path for a given interest rate path. wpath, divpath, etc. are guesses for the transition path
"""

function solve_for_transition(Rpath::Array{Float64,1},wpath::Array{Float64,1},div_path::Array{Float64,1},Spath::Array{Float64,1},SS::steady_state,p::params)

#incomplete

end



"""
    simulate_forward()

Simulates a distribution of households and computes total consumption (C), labor supply (L), asset position (B) and wealth distribution (D)
"""

function simulate_forward(D0::Array{Float64,1},cpol_path::Array{Float64,2},Rpath::Array{Float64,1},
                        wpath::Array{Float64,1},div_path::Array{Float64,1},τ_path::Array{Float64,1},p::params)

 #back out number of periods
 T = length(Rpath)

 #pre-allocate some arrays
 Cpath = Array{Float64,1}(undef,T-1);  Lpath = Array{Float64,1}(undef,T-1)
 Bpath = Array{Float64,1}(undef,T-1);  Dpath = Array{Float64,2}(undef,p.nz*p.nk,T-1)

 Dpath[:,1] = D0

 for t = 2:T-1
 Cpath[t], Lpath[t], Bpath[t], Dpath[:,t] = simulate_step(Dpath[:,t-1],reshape_c(cpol_path[:,t],p)
                                            ,Rpath[t],wpath[t],τ_path[t],div_path[t],p)
 end
 

 return Cpath[2:end],Lpath[2:end]

    
end

"""
    simulate_step()

Conducts forward simulation for one period. Helper function to simulate_forward, equivalent to simulatestep() in MNS code.

"""
function simulate_step(D::Array{Float64,1},c_pol::Array{Float64,2},R::Float64,w::Float64,τ::Float64,div::Float64,p::params)

    @unpack k_grid = p

 #household consumption and labor supply
 C,L = aggregate_C_L(D,c_pol,R,w,τ,div,p)

 Pi = forwardmat(c_pol,R,w,τ,div,p)

 #distribution in next period
 Dprime = Pi*D

 #aggregate assets
 Assets = dot(Dprime,repeat(k_grid,3,1))

 #test for validity of distribution 
 @assert (abs(sum(D) - 1.0) < 1e-6)

 return C,L,Dprime,Assets
    
end