#features functions to calibrate and compute the steady state

using Interpolations, SparseArrays, LinearAlgebra, NLsolve

"""
   get_steady_state()

Finds the steady state for given interest rate by iterating over β, using a simple bisection type updating algorithm.
I typically used Y = 0.6 as guess for Y and beta_Range = [0.95,0.99], which works fine. Y_guess not particularly important.
"""

#find steady state equilibrium using different method
function get_steady_state(p::params, Y_guess::Float64, beta_Range::Array{Float64,1};
                           tolβ::Float64 = 1e-4, tolY::Float64=1e-6  )

   @unpack μ,Rbar,B,Γ,tax_weights,b_grid,k_grid = p

   #outer loop: bisection to find steady state interest rate
   iterβ = 1; distβ = 1.0; distY = 1.0
   β = beta_Range[2]; Ys = Y_guess 
   βmax = beta_Range[2]; βmin = beta_Range[1]

   #generate objects so that they exist outside of while scope
   agg_assets = 0.0; D = 0.0; C=0.0; L = 0.0;
   w = 0.0; τ=0.0; div = 0.0; 

   #initial guess for consumption policy
   c_pol = 0.3 .+ 0.1*p.b_grid; c_pol = repeat(c_pol,1,3)

   while (iterβ<200) & ( (distβ > tolβ) | (distY > tolY) )
      
      #distY = 1.0; iterY = 1
      #while (iterY < 5) & (distY>tolY)
         #inner loop will initially not converge

         R = Rbar; w = 1/μ; τ = B*Ys*(1-1/R)/(Γ'*tax_weights); div = Ys*(1-w)
         
         #get policy functions
         c_pol = EGM_SS(c_pol, β, Ys, p)

         #get wealth-income transition matrix
         Pi = forwardmat(c_pol,R,w,τ,div,p)

         #get Steady State wealth distribution
         D = inv_dist(Pi)
    
         #compute aggregate assets
         agg_assets = dot(D,repeat(k_grid,3,1))

         #get aggregate comsumption and labor (L = steady state output)
         C, L = aggregate_C_L(D,c_pol,R,w,τ,div,p)

         distY = abs(C-L)
         #Ys = 0.9*Ys+0.1*C #heuristic update rule
         Ys = C

       #  iterY = iterY+1
         #println("Y-iteration: ",iterY," Distance: ",distY," Current Y:",Ys)
      #end
    

    distβ = abs(agg_assets - B*Ys)
    println("β-iteration: ",iterβ," Dist. Assets: ",distβ," Dist. Y: ",distY, " Current β: ", β)
    println(" ")

    #update β according to bisection rule
    if β == beta_Range[2]
      @assert (agg_assets > B*Ys)
      β = beta_Range[1]
    elseif β == beta_Range[1]
      @assert (agg_assets < B*Ys)
      βprev = β[1]
      β = mean(beta_Range)
    else 
         if (agg_assets > B*Ys)
         βmax = β[1]
         elseif (agg_assets < B*Ys)
         βmin = β[1]
         end
         β = mean([βmax,βmin])
    end
      
    iterβ = iterβ+1
   end

   SS_structure = steady_state(c_pol,D,Ys,C,L,agg_assets,w,τ,div,Rbar)

 return β, Ys, SS_structure
end


"""
    check_steady_state(β::Float64,Y::Float64,p::params)

For given β,Y and Parameters, this function computes consumption and labor suppy 
as implied by the household decisions and returns the distance between output 
implied by labor and output implied by consumption as well as the distance between household
asset choices and the asset target. It can be used in an Root-Finding algorithm to find the β
corresponding to the asset and interest rate target.
"""
function check_steady_state(beta::Float64,Y::Float64,p::params;return_distance::Bool=false)

    @unpack μ,Rbar,B,Γ,tax_weights,b_grid,k_grid = p

    β = beta/100.0
    
    #steady state prices:
    R = Rbar; w = 1/μ; τ = B*Y*(1-1/R)/(Γ'*tax_weights); div = Y*(1-w)

    #initial guess for consumption policy 
    c_guess = 0.3 .+ 0.1*p.b_grid; c_guess = repeat(c_guess,1,3)
    
    #find HH policy function
    c_policies = EGM_SS(c_guess, β, Y, p)
   
    #get wealth-income transition matrix
    Pi = forwardmat(c_policies,R,w,τ,div,p)

    #get Steady State wealth distribution
    D = inv_dist(Pi)
    
    #compute aggregate assets
    agg_assets = dot(D,repeat(k_grid,3,1))

    #get aggregate comsumption and labor (L = steady state output)
    C, L = aggregate_C_L(D,c_policies,R,w,τ,div,p)

    #compute distance vector 
    SS_dist = Array{Float64,1}(undef,2)
    SS_dist[1] = agg_assets/C - B 
    SS_dist[2] = C - L
   
    if return_distance == true
    println("Current C: ",C," Current K: ",agg_assets)
    println("Current Dist[1]: ",SS_dist[1]," Current Dist[2]: ",SS_dist[2])
    println(" ")
    return SS_dist
    else
    return SS_structure = steady_state(c_policies,D,Y,C,L,agg_assets,w,τ,div,R)
    end

end


"""
    EGM_SS(β::Float64,Y::Float64,p::params)

Computes household steady state (SS) savings and labor supply using the endogenous grid method (EGM).
"""
function  EGM_SS(c_guess::Array{Float64,2},β::Float64,Y::Float64,p::params;
                 T::Int=30, maxit::Int = 500, tol::Float64 = 1e-7)

   @unpack μ,Rbar,B,Γ,tax_weights = p

   #steady state prices:
   R = Rbar; w = 1/μ; τ = B*Y*(1-1/R)/(Γ'*tax_weights); D = Y*(1-w)
    
    iter = 0; dist = 1 #initialize iteration
    while (iter < maxit) & (dist > tol)  

    c_policies = solveback(reshape_c(c_guess,p),w*ones(T),R*ones(T),τ*ones(T),
    D*ones(T),β,p)

    dist = maximum(abs.(c_policies[:,1].-c_policies[:,2])) #measure maximum distance

    c_guess = reshape_c(c_policies[:,1],p) #update policy rules

    iter = iter + 1 #count iterations
    end
    #println("Steady State Analytics")
    println("No. of iterations: ",iter,"  Distance: ",dist)
 
   return c_guess #return steady state policy functions

end


"""
    solveback()

 Solves backwards from a period (e.g. steady state) in which decision rules are known
 and given price/dividend paths. Relies on function EGM.
 (at the moment, this still omits β-heterogeneity).
 In w_path etc, the final values must correspond to the period to solve back from
"""

function solveback(c_final::Array{Float64,1},w_path::Array{Float64,1},R_path::Array{Float64,1},τ_path::Array{Float64},
                    div_path::Array{Float64,1},β::Float64,p::params) 

  @unpack nb, nz = p

  #calculate total amount of periods
  T = size(w_path)[1]

  #pre-allocate some arrays for consumption paths
  c_path = zeros(size(repeat(c_final,1,T)))
  temp   = Array{Float64,2}(undef,nb,nz)

  #set final values
  c_path[:,end] .= c_final 

  #solving back
  for t = T-1:-1:1

  temp[:,:] .=  EGM(reshape_c(c_path[:,t+1],p),β, R_path[t:t+1], w_path[t:t+1],
               τ_path[t:t+1], div_path[t:t+1],p)

   c_path[:,t] .= reshape_c(temp,p)

  end
    
  return c_path
    
end


"""
    margU(c::Float64,par::params)

 Helper function, computes marginal utility of consumption or second derivative of utility function (if order = 2)
"""
function margU(c::Float64,par::params,order::Int=1)
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
w::Float64,τ::Float64,div::Float64,inc_idx::Int,p::params)

 @unpack b_grid,ψ,z,tax_weights = p

 itp = LinearInterpolation(b_grid,cons[:,inc_idx],extrapolation_bc = Line()) #create interpolant

 c = itp.(xthis) #interpolate


 n = (margU.(c,(p,)).*w.*z[inc_idx]).^(1/ψ)
 bp = R*(xthis.+w.*n.*z[inc_idx] .- τ*tax_weights[inc_idx].+ div .- c )

 return c,n,bp
end

"""
    EGM(β::Float64,...)

Conducts one iteration on consumption and labor supply using the endogenous grid method (EGM).
Inputs need to be 2x1 arrays of wage, interest etc in current and previous period
""" 
function EGM(c_next::Array{Float64,2},β::Float64, Rs::Array{Float64,1}, ws::Array{Float64,1},
                τs::Array{Float64,1}, div::Array{Float64,1},p::params)

 @unpack nb, nz, b_grid, γ, z, Πz, ψ, tax_weights = p 

 #asset grid to use
 xthis = vcat(0.0,b_grid)
 #size asset grid
 nx = nb + 1

 #pre-allocate arrays for marginal utilities
 MU = Array{Float64,2}(undef,nx,nz)

 #calculate marginal utilities
 for i = 1:nz

    cthis, = get_cnbp(xthis, c_next,Rs[2],ws[2],τs[2],div[2],i,p)

    @assert (sum(cthis .> 0.0) .== nx) #check whether all consumption values are positve

    MU[:,i] .= margU.(cthis,(p,))
 end
 
 #compute expected marginal utilities for the different income types
 MUexp = Array{Float64,2}(undef,nx,nz)  #pre-allocate
 for i = 1:nz
 MUexp[:,i] .= sum(MU.*Πz[i,:]',dims=2)[:,1]
 end

 Cprev = (β*Rs[1]*MUexp).^(-(1/γ))

 @assert all(x -> x>0, Cprev) #check whether all consumption values are positve

 labor = (margU.(Cprev,(p,)).*ws[1].*z').^(1/ψ)

 bprev = xthis./Rs[1] .+ Cprev .- ws[1].*labor.*z' .+ τs[1].*tax_weights' .- div[1]
 
 #generate c_new
 Cnew = Array{Float64,2}(undef,nb,nz)
 
 for i = 1:nz
    ind_constrained = (b_grid .<= bprev[1,i])
    
    if any(ind_constrained)
        Cnew[ind_constrained,i].= egm_solve_constrained(b_grid[ind_constrained],i,ws[1],τs[1],div[1],Rs[1],p)
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
function egm_solve_constrained(bs::Array{Float64,1},ip::Int,w::Float64,τ::Float64,div::Float64,R::Float64,p::params)

 #initial guess for labor supply
 n = 0.6.*ones(size(bs))
 c = bs .+ n.*w*p.z[ip] .+ div .+ (1-1/R)*p.a_min .- τ*p.tax_weights[ip]

 iter = 0; dist = 1.0
 while (iter < 1000) & (dist > 1e-7)  

    f = margU.(c,(p,)).*w.*p.z[ip] .- n.^p.ψ
    J = margU.(c,(p,),2).*(w*p.z[ip])^2 .- p.ψ*n.^(p.ψ-1.0)
    
    #update 
    n = n .- f./J
    c = bs .+ n.*w*p.z[ip] .+ div .+ (1-1/R)*p.a_min .- τ*p.tax_weights[ip]

    dist = maximum(abs.(f))
    iter = iter+1
 end

 #check for convergence
 if iter == 1000
    error("egm_solve_constrained did not converge")
 end

 c_cons = c

 return c_cons 

end


"""
   forwardmat(c_opt::Array{Float64,2},R::Float64,w::Float64,τ::Float64,div::Float64)

Generates a transition matrix for the aggregate wealth distribution. Uses the function lineartrans() as 
in the original code.
"""

function forwardmat(c_opt::Array{Float64,2},R::Float64,w::Float64,τ::Float64,div::Float64, par::params)

   @unpack nk,nz,k_grid,Πz = par


   #pre-allocate some Arrays:
   IR = Array{Int64,1}(undef,2*nk*nz^2); IC = Array{Int64,1}(undef,2*nk*nz^2)
   VV =  Array{Float64,1}(undef,2*nk*nz^2)
   
   idx1 = 1 #variable used for indexing below
   for i = 1:nz
   
    bp = get_cnbp(k_grid,c_opt,R,w,τ,div,i,par)[3] #saving choices for income type

    From, To, Probs = lineartrans(bp,par) #get transitions

    #aux. variables used for indexing below
    offsi = (i-1)*nk 

      for j = 1:nz

         ppj = Πz[i,j]

         #aux. variable used for indexing below
         offsj = (j-1)*nk; idx2 = idx1 + 2*nk - 1
         
         #this will generate the row indices of the sparse transition matrix
         #this indicates the position the agents come from
         IR[idx1:idx2] = offsi .+ From 
         #this will generate the column indices of the sparse transition matrix
         #this where they go to
         IC[idx1:idx2] = offsj .+ To
         #this will be the entries of the sparse transition matrix
         #and this the probabilities where they will go to
         VV[idx1:idx2] = ppj*Probs

             
         idx1 = idx1 + 2*nk #update index variable
      end
   

   end

 return Pi = sparse(IR,IC,VV,nz*nk,nz*nk) #return transition matrix
end



"""
   lineartrans(bp::Array{Float64,1},par::params)

Calculates the new positions and transition probabilities on the asset grid for given savings choices.
Used the to conduct a non-stochastic simulation a la Young (2010, JEDC).
"""
function lineartrans(bp::Array{Float64,1},par::params)
 @unpack nk, k_grid = par
 
 #remove values higher than maximum on asset grid
 k0 = minimum(hcat(bp,maximum(k_grid)*ones(size(bp))),dims = 2)[:]

 #make sure all values are above 0
 k0 = maximum(hcat(k0,zeros(size(k0))),dims = 2)[:]

 #find next closest value on grid that is lower for each element of k0
 iPos = find_closest_lower.(k0,(k_grid,))
 #make sure no value in iPos exceed nk - 1
 iPos = minimum(hcat(iPos,(nk-1)*repeat([1],nk)),dims = 2)[:] 
 
 
 pHigh = (k0 .- k_grid[iPos])./(k_grid[iPos.+1]-k_grid[iPos])

 From = repeat(collect(1:nk),2,1)
 To   = vcat(iPos,iPos .+ 1)
 Probs= vcat(1 .- pHigh,pHigh)

 return From, To, Probs

end


"""
   find_closest_lower(x::Float64,grid::Array{Float64,1})

Helper function for lineartrans
"""
function find_closest_lower(x::Float64,grid::Array{Float64,1}=k_grid)

   closest_index = argmin(abs.(x.-grid))
   if grid[closest_index] <= x
      return closest_index
   else
      return closest_index - 1
   end
end

"""
   inv_dist

Find invariant distribution by iterating over transition matrix
uses insights from https://discourse.julialang.org/t/stationary-distribution-with-sparse-transition-matrix/40301
"""
function inv_dist(Π::SparseMatrixCSC)

   x = [1; (I - Π'[2:end,2:end]) \ Vector(Π'[2:end,1])]
   return dist = x./sum(x) 

   #uncomment below for conventional method - takes longer
   
   #n = size(Π)[1]
   #mc = MarkovChain(Array(Π),collect(1:n))
   #return stationary_distributions(mc)[1]

end


"""
   aggregate_C_L(D::Array{Float64,1},c_policies::Array{Float64,1},R::Float64,w::Float64,τ::Float64,div::Float64)

Computes aggregate consumption and labor supply of the household sector. Equivalent to expect_C in the MNS code.
"""
function aggregate_C_L(D::Array{Float64,1},c_policies::Array{Float64,2},R::Float64,w::Float64,τ::Float64,div::Float64,p::params)
  
   @unpack nz,nk,k_grid,z = p

   #loop over income states
   L = 0.0; C = 0.0; #initialize
   for i = 1:nz

      #for simple indexing
      idx = (i - 1)*nk

      #get consumption and labor supply for income/wealth bin
      c,n, = get_cnbp(k_grid,c_policies,R,w,τ,div,i,p)

      #sum up (using dot product)
      C = C + dot(D[idx+1:idx+nk],c)
      L = L + dot(D[idx+1:idx+nk],n*z[i])

   end

   return C,L
end


"""
   get_steady_state(p::params,β_guess::Float64,Y_guess::Float64)

Solve for β and steady state output so target bond level is reached at the target interest rate (both are specified in p)
DOES NOT WORK YET!

function get_steady_state(p::params,β_guess::Float64,Y_guess::Float64)

   #@unpack  = p
   
   #define objective for nonlinear solver
   function obj!(Dist,x)
      Dist = check_steady_state(x[1],x[2],p;return_distance=true)
   end

   SR = nlsolve(obj!,[β_guess*100,Y_guess]  ) #get solver results

   #return steady state structure and the β value we solved for 
   return SR

end
"""