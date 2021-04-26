#features functions to calibrate and compute the steady state

using Interpolations

"""
    check_steady_state(β::Float64,Y::Float64,p::params)

For given β,Y and Parameters, this function computes consumption and labor suppy 
as implied by the household decisions and returns the distance between output 
implied by labor and output implied by consumption as well as the distance between household
asset choices and the asset target. It can be used in an Root-Finding algorithm to find the β
corresponding to the asset and interest rate target.
"""
function check_steady_state(β::Float64,Y::Float64,p::params)

    @unpack μ,Rbar,B,Γ,tax_weights = p

    #steady state prices:
    R = Rbar; w = 1/μ; τ = B*Y*(1-1/R)*(Γ'*tax_weights); D = Y*(1-w)

    #not complete yet

end


"""
    EGM_SS(β::Float64,Y::Float64,p::params)

Computes household steady state (SS) savings and labor supply using the endogenous grid method (EGM).
"""
function  EGM_SS(c_guess::Array{Float64},β::Float64,Y::Float64,p::params;
                 T::Int=30, maxit::Int = 500, tol::Float64 = 1e-7)
    
    iter = 0; dist = 1 #initialize iteration
    while (iter < maxit) & (dist > tol)  

    c_policies = 0 #solveback here
    
    dist = maximum(abs.(c_policies[:,1].-c_policies[:,2]))

    c_guess = c_policies[:,1] #update policy rules
    end
 
    return c_policies[:,1] #return steady state policy functions

end


"""
    solveback()

 Solves backwards from a period (e.g. steady state) in which decision rules are known
 and given price/dividend paths. Relies on function EGM.
 (at the moment, this still omits β-heterogeneity).
 In w_path etc, the final values must correspond to the period to solve back from
"""

function solveback(c_final::Array{Float64,1},w_path::Array{Float64,1},τ_path::Array{Float64},
                    div_path::Array{Float64,1},β::Float64,p::params) 

  #calculate total amount of periods
  T = size(w_path)[1]

  #pre-allocate some arrays for consumption paths
  c_path = zeros(size(repeat(c_final,1,T)))
  #set final values
  c_path[:,end] = c_final 

  #solving back
  for t = T-1:-1:1

    c_path[:,t] = 0 # here EGM function

  end
    

    
end


"""
    margU(c::Float64,par::params=p)

 Helper function, computes marginal utility of consumption or second derivative of utility function (if order = 2)
"""
function margU(c::Float64,order::Int=1,par::params=p)
 if order == 1
    return mu = c^(-par.γ)
 elseif order == 2
    return -par.γ*c^(-par.γ-1)
 else
    error("Order must be either 1 or 2!")
 end
end


"""
    get_cnbp(xthis::Array{Float64},c::Array{Float64},R::Float64,
            w::Float64,τ::Float64,p::params,inc_idx::Int)

Helper function to retrieve some values

"""
function get_cnbp(xthis::Array{Float64},cons::Array{Float64,2},R::Float64,
w::Float64,τ::Float64,div::Float64,inc_idx::Int,p::params=p)

 @unpack b_grid,ψ,z,tax_weights = p

 itp = LinearInterpolation(b_grid,cons[:,inc_idx]) #create interpolant

 c = itp.(xthis) #interpolate


 n = (margU.(c).*w.*z[inc_idx]).^(1/ψ)
 bp = R*(xthis.+wage.*n.*z[inc_idx] .- τ*tax_weights[inc_index].+ div .- c )

 return c,n,bp
end

"""
    EGM(β::Float64,...)

Conducts one iteration on consumption and labor supply using the endogenous grid method (EGM).
Inputs need to be 2x1 arrays of wage, interest etc in current and previous period
""" 
function EGM(c_next::Array{Float64,2},β::Float64, Rs::Array{Float64,1}, ws::Array{Float64,1},
                τs::Array{Float64,1}, div::Array{Float64,1},p::params)

 @unpack nb, nz, b_grid, γ, z = p 

 #asset grid to use
 xthis = vcat(0.0,b_grid)
 #size asset grid
 nx = nb + 1

 #pre-allocate arrays for marginal utilities
 MU = Array{Float64,2}(undef,nx,nz)

 #calculate marginal utilities
 for i = 1:nz

    cthis, = get_cnbp(xthis, c_next,R[2],w[2],τs[2],div[2],i)

    @assert (sum(cthis .> 0.0) .== nz) #check whether all consumption values are positve

    MU[:,i] = margU.(cthis)
 end
 
 #compute expected marginal utilities for the different income types
 MUexp = Array{Float64,2}(undef,nx,nz)  #pre-allocate
 for i = 1:nz
 MUexp[:,i] = MU[:,i].*Πz[i,:]
 end

 Cprev = (β*Rs[1]*margU.(MUexp)).^(-(1/γ))

 @assert all(x -> x>0, Cprev) #check whether all consumption values are positve

 labor = (margU.(Cprev).*ws[1].*z').^(1/ψ)

 bprev = xthis./Rs[1] .+ Cprev .- ws[1].*labor.*z' + τs[1]*tax_weights' - div[1]

 #generate c_new
 Cnew = Array{Float64,2}(undef,nx,nz)

 for i = 1:nz
    ind_constrained = (b_grid .<= bprev[1,i])

    if any(ind_constrained)
        Cnew[ind_constrained,i].= egm_solve_constrained(b_grid[ind_constrained],i,ws[1],τs[1],div[1],Rs[1])
    end
 itp = LinearInterpolation(bprev[:,i],Cprev[:,i],extrapolation_bc = Line())

 Cnew[.!ind_constrained,i] .= itp(b_grid[.!ind_constrained]) 
 end

 return Cnew #return results
end


"""
    egm_solve_constrained()

 Backs out consumption level of constrained household.
"""
function egm_solve_constrained(bs::Array{Float64,1},ip::Int,w::Float64,τ::Float64,div::Float64,R::Float64,p::params=p)

 #initial guess for labor supply
 n = 0.6.*ones(size(bs))
 c = b .+ n.*w*p.z[ip] .+ div .+ (1-1/R)*p.a_min

 iter = 0; dist = 1.0
 while (iter < 1000) & (dist > 1e-7)  

    f = margU.(c).*w.*p.z[ip] .- n.^p.ψ
    J = margU(c,2).*(w*p.z[ip])^2 .- p.ψ*n.^(p.ψ-1.0)
    
    #update 
    n = n .- f./J
    c = b .+ n.*w*p.z[ip] .+ div .+ (1-1/R)*p.a_min

    dist = abs(f)
    iter = iter+1
 end

 #check for convergence
 if iter == 1000
    error("egm_solve_constrained did not converge")
 end

end
