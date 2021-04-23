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
    get_cnbp(xthis::Array{Float64},c::Array{Float64},Rnext::Float64,
            wnext::Float64,τnext::Float64)

Helper function to retrieve some values

"""
function get_cnbp(xthis::Array{Float64},c::Array{Float64},Rnext::Float64,
wnext::Float64,τnext::Float64,p::params)

return c,n,bp
end

"""
    EGM(β::Float64,...)

Conducts one iteration on consumption and labor supply using the endogenous grid method (EGM).
Inputs need to be 2x1 arrays of wage, interest etc in current and previous period
""" 
function EGM(c_next::Array{Float64,2},β::Float64, Rs::Array{Float64,1}, ws::Array{Float64,1},
                τs::Array{Float64,1}, div::Array{Float64,1},p::params)

 @unpack nb, nz, b_grid = p 

 #asset grid to use
 xthis = vcat(0.0,b_grid)
 #size asset grid
 nx = nb + 1

 #pre-allocate arrays for marginal utilities
 MU = Array{Float64,2}(undef,nx,nz)


 #incomplete




end