

using Parameters, Roots, QuantEcon, Statistics


# File that provides structures used throughout the code as well as functions to fill them

"""
    A structure that stores the parameter vector of the model, including all grids.
    It can be automatically generayed using function `set_parameters()`.
    It contains the following elements: 
    ..* The discount factor `β` 
    ..* The CRRA risk aversion parameter `γ` 
    ..* Yhe inverse Frisch elasticity `ψ`
    ..* The labor disutility scaling parameter `ψ1` (only used in complete markets case)
    ..* The target asset level `B`
    ..* The markup/CES production function parameter `μ`
    ..* The Calvo price re-set probability `θ`
    ..* The target interest rate `Rbar`
    ..* The number of grid points on the labor productivity grid `nz`
    ..* The labor productivity grid `z`
    ..* The labor productivity state Markov transition matrix `Πz`
    ..* The invariant distribution over labor productivty states `Γ`
    ..* The tax weights `tax_weights`
    ..* The borrowing constraint `a_min`
    ..* The maximum value on the asset grid `a_max`
    ..* The number of points on the asset grid used for computing the wealth distribution `nk` 
    ..* The number of points on the asset grid used for computing HH policy functions `nb` 
    ..* The asset grid used for computing the wealth distribution `k_grid`
    ..* The asset grid used for computing the HH policy functions `b_grid`

"""
 @with_kw mutable struct  params
    β::Float64 = 0.986    #discount factor
    γ::Float64 = 2.0      #Risk aversion
    ψ::Float64  = 2.0     #inverse Frisch elasticity
    ψ1::Float64 = 1.0     #scaling labor disutility (only used in CompMkts part)
    B::Float64 = 5.5      #Supply of Assets: approx. 1.4 times steady state GDP
    μ::Float64 = 1.2      #Markup
    θ::Float64 = 0.15     #Probability that Calvo fairy visits
    Rbar::Float64 = 1.005 #Target quarterly interest rate
    nz::Int64               #Number of income grid points
    z::Array{Float64,1}   #Possible income states
    Πz::Array{Float64,2}           #Markov transition matrix
    Γ::Array{Float64,1}            #invariant distribution over income states
    tax_weights::Array{Float64,1}  #tax weights
    a_min::Float64        #minimum asset holdings in grid
    a_max::Float64        #maximum asset holdings in grid
    nk::Int64               #number of grid points used for quadrature grid
    nb::Int64               #number of consumption polynomial knot points
    k_grid::Array{Float64,1}       #quadrature grid for wealth distribution
    b_grid::Array{Float64,1}       #knot points for consumption polynomial
end


""" 
    A structure that saves the steady state of the incomplete markets model.
    It contains the household policy functions (`c_polcicies`), the invariant distribution of wealth (`D`),
    output (`Y`), consumption (`C`), hours worked (`L`), aggregate assets (`K`), the wage (`w`), the tax rate (`τ`),
    the dividebd (`div`) and the interest rate (`Rbar`). Note that `C`, `L` and `Y` should be equal (approximately).
    
"""
mutable struct steady_state
 c_policies::Array{Float64}  #policy functions
 D::Array{Float64,1}         #invariant distribution over income/wealth types
 Y::Float64                  #Aggregate output
 C::Float64                  #Aggregate consumption
 L::Float64                  #Aggregate labor supply
 K::Float64                  #aggregate assets 
 w::Float64                  #wage
 τ::Float64                  #taxes
 div::Float64                #dividends
 R::Float64                  #interest rate
end


"""
    A structure that stores the incomplete markets model transition path, including the time paths of price dispersion (`S`), the wage (`w`),
    inflation (`pΠ`), output (`Y`), the interest rate (`R`), the tax rate (`τ`) and the dividend (`div`).
"""
mutable struct transition_full
S::Array{Float64,1}         #price dispersion
w::Array{Float64,1}         #wages
pΠ::Array{Float64,1}        #reset inflation
Y::Array{Float64,1}         #output
R::Array{Float64,1}         #interest rates
τ::Array{Float64,1}         #tax rates
div::Array{Float64,1}       #dividends
end   


"""
    A structure that stores the complete markets transition path, including the time paths of price dispersion (`S`), the wage (`w`),
    inflation (`pΠ`), output (`Y`), the interest rate (`R`),  and the dividend (`div`).
"""
struct transition_CompMkts
S::Array{Float64,1}         #price dispersion
w::Array{Float64,1}         #wages
pΠ::Array{Float64,1}        #reset inflation
Y::Array{Float64,1}         #output
R::Array{Float64,1}         #interest rates
div::Array{Float64,1}       #dividends
end


"""
    logspaceshift(xa::Float64,xb::Float64,n::Int,x2::Float64,n_at_x2::Float64)

A function involved in replicating a log-spaced asset grids exactly as in the MNS code. Several methdos available.
"""
function logspaceshift(xa::Float64,xb::Float64,n::Int,x2::Float64,n_at_x2::Float64)
    frac = n_at_x2 / n

    f(x) = log(x2+x)-log(xa+x) - frac*(log(xb+x)-log(xa+x)) #auxiliary function
    xshift = find_zero(f, 0.5) 

    grid = exp.(collect(range(log(xa+xshift),log(xb+xshift),length=n))).-xshift
    grid[1] = xa ; grid[end] = xb

    return grid, xshift
 
end

"""
    logspaceshift(xa::Float64,xb::Float64,n::Int,x2::Float64)

A function involved in replicating a log-spaced asset grids exactly as in the MNS code. Several methdos available.
This is one of them.
"""
function logspaceshift(xa::Float64,xb::Float64,n::Int,x2::Float64)

    xshift = x2
    grid = exp.(collect(range(log(xa+xshift),log(xb+xshift),length=n))).-xshift
    grid[1] = xa ; grid[end] = xb

    return grid
end


"""
    logspaceshift(xa::Float64,xb::Float64,n::Int)

A function involved in replicating a log-spaced asset grids exactly as in the MNS code. Several methdos available.
This is one of them.
"""
function logspaceshift(xa::Float64,xb::Float64,n::Int)

    xshift = 0.0
    grid = exp.(collect(range(log(xa+xshift),log(xb+xshift),length=n))).-xshift
    grid[1] = xa ; grid[end] = xb

    return grid, xshift
end


"""
    makeknotd(kmin::Float64,kmax::Float64,n::Int,logshift::Float64)

A Function that produces an asset grid using the function `logspaceshift`. Several methods available. This is one of them.
"""
function makeknotd(kmin::Float64,kmax::Float64, n::Int, logshift::Float64)

    knots = logspaceshift(kmin,kmax,n,logshift)
    return knots
    
end 

"""
    makeknotd(kmin::Float64,kmax::Float64,n::Int,logshift::Float64)

A Function that produces an asset grid using the function `logspaceshift`. Several methods available. This os one of them.
"""
function makeknotd(kmin::Float64,kmax::Float64, n::Int)

    knots,logshift = logspaceshift(kmin,kmax,n,1.0,n/4)
    return knots, logshift
end 


"""
    set_parameters(;β::Float64 = 0.986,    # Discount factor
                    γ::Float64 = 2.0,      # Risk aversion
                    ψ::Float64 = 2.0,      # Inverse Frisch elasticity
                    ψ1::Float64 = 1.0,     # Labor disutility scaling (only in complete markets)
                    B::Float64 = 5.5,      # Supply of Assets: 1.4 times annual GDP = 5.6 times quarterly GDP
                    μ::Float64 = 1.2,      # Mark-up
                    θ::Float64 = 0.15,     # Probability that Calvo fairy visits
                    Rbar::Float64 = 1.005, # Target quarterly interest rate
                    nz::Int64 = 3,         # Number of income grid states
                    ρ::Float64 = 0.96566,            # Persistence parameter income process
                    σ::Float64 = 0.01695^0.5,        # Std of income process random comp
                    tax_weights::Array{Float64,1} = [0.0,0.0,1.0],   # Tax weights
                    a_min::Float64  = 0.0,           # Borrowing constraint
                    a_max::Float64 =75.0,            # Maximum asset holdings in grid
                    x_min::Float64 = 0.001,          # Parameter to generate the consumption knots
                    nk::Int64 = 1000,                # Number of grid points used for quadrature grid
                    nb::Int64 = 200)

Function to construct the parameter structure. The default parameters correspond to the baseline calibration from the MNS paper.
The function applies the Rouwenhorst method (using its `QuantEcon` implementation) to discretize the AR(1) process for (log) labor productivity
with persistence ρ and innovation variance σ.    
Note that while `\beta` is set automatically, it will typically be replaced after the steady state is computed, since
β is numerically computed to match a target interest rate and aggregate asset level.

To generate the asset grids, the functions `makeknotd` and `logspaceshift` are used.
"""
function set_parameters( ; β::Float64 = 0.986,    #discount factor
                           γ::Float64 = 2.0,       #Risk aversion
                           ψ::Float64 = 2.0,      #inverse Frisch elasticity
                           ψ1::Float64 = 1.0,     #labor disutility scaling (only in complete markets)
                           B::Float64 = 5.5,      #Supply of Assets: 1.4 times annual GDP = 5.6 times quarterly GDP
                           μ::Float64 = 1.2,      #Markup
                           θ::Float64 = 0.15,     #Probability that Calvo fairy visits
                           Rbar::Float64 = 1.005, #Target quarterly interest rate
                           nz::Int64 = 3,           #Number of income grid states
                           ρ::Float64 = 0.96566,            #persistence parameter income process
                           σ::Float64 = 0.01695^0.5,          #Std of income process random comp.
                           tax_weights::Array{Float64,1} = [0.0,0.0,1.0],   #tax weights
                           a_min::Float64  = 0.0,            #borrowing constraint
                           a_max::Float64 =75.0,             #maximum asset holdings in grid
                           x_min::Float64 = 0.001,          #parameter to generate the consumption knots
                           nk::Int64 = 1000,                  #number of grid points used for quadrature grid
                           nb::Int64 = 200 )                   #number of consumption polynomial knot points

 #check whether tax_weights have consistent size 
 @assert (size(tax_weights,1) == nz) "Length of tax weights must equal number of income grid points"

 #obtain discrete income process using Rouwenhorst (1995) method
 mc = rouwenhorst(nz,ρ,σ)
 #get stationary distribution
 Γ = stationary_distributions(mc)[1]
 #extract transition matrix and income grid points
 z = exp.(collect(mc.state_values)); Πz = mc.p
 
 #get grids in the same way as in paper
 k_grid, logshift = makeknotd(a_min,a_max,nk)

 b_grid = logspaceshift(x_min,a_max,nb,logshift)

 return params(β,γ,ψ,ψ1,B,μ,θ,Rbar,nz,z,Πz,Γ,tax_weights.*z,a_min,a_max,nk,nb,k_grid,b_grid)

end


"""
    reshape_c(c::Array{Float64},p::params)

Helper functions that gets the HH consumption policy functions into a different form. If they are supplied as a `nb*nz` matrix, they 
will be turned into a vector of length `nb+nz` and vice versa. This is the method doing the latter.
"""
function reshape_c(c::Array{Float64,1},par::params)
 @unpack nb,nz = par

 c_wide = reshape(c,(nb,nz))

 return c_wide
 
end

"""
    reshape_c(c::Array{Float64,2},par::params)

Helper functions that gets the HH consumption policy functions into a different form. If they are supplied as a `nb*nz` matrix, they 
will be turned into a vector of length `nb+nz` and vice versa. This is the method doing the former.
"""
function reshape_c(c::Array{Float64,2},par::params)
 @unpack nb,nz = par

 c_long = reshape(c,(nb*nz,1))

 return c_long[:,1]

end