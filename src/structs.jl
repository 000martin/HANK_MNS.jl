

using Parameters, Roots, QuantEcon, Statistics


#File that provides structures used throughout the code as well as functions to fill them

"""
A structure that stores the parameter vector of the model
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
    nz::Int               #Number of income grid points
    z::Array{Float64,1}   #Possible income states
    Πz::Array{Float64,2}           #Markov transition matrix
    Γ::Array{Float64,1}            #invariant distribution over income states
    tax_weights::Array{Float64,1}  #tax weights
    a_min::Float64        #minimum asset holdings in grid
    a_max::Float64        #maximum asset holdings in grid
    nk::Int               #number of grid points used for quadrature grid
    nb::Int               #number of consumption polynomial knot points
    k_grid::Array{Float64,1}       #quadrature grid for wealth distribution
    b_grid::Array{Float64,1}       #knot points for consumption polynomial
end

""" 
A structure that saves the steady state of the model.
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
A structure that stores the full model transition
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
A structure that stores the complete markets transition
"""
struct transition_complete_markets
S::Array{Float64,1}         #price dispersion
w::Array{Float64,1}         #wages
pΠ::Array{Float64,1}        #reset inflation
Y::Array{Float64,1}         #output
R::Array{Float64,1}         #interest rates
div::Array{Float64,1}       #dividends
end



"""
    logspaceshift(xa::Float64,xb::Float64,n::Int,x2::Float64,n_at_x2::Float64)

A function involved in making the grids. Several methods available.
"""

function logspaceshift(xa::Float64,xb::Float64,n::Int,x2::Float64,n_at_x2::Float64)
    frac = n_at_x2 / n

    f(x) = log(x2+x)-log(xa+x) - frac*(log(xb+x)-log(xa+x)) #auxiliary function
    xshift = find_zero(f, 0.5) 

    grid = exp.(collect(range(log(xa+xshift),log(xb+xshift),length=n))).-xshift
    grid[1] = xa ; grid[end] = xb

    return grid, xshift
 
end

function logspaceshift(xa::Float64,xb::Float64,n::Int,x2::Float64)

    xshift = x2
    grid = exp.(collect(range(log(xa+xshift),log(xb+xshift),length=n))).-xshift
    grid[1] = xa ; grid[end] = xb

    return grid
end

function logspaceshift(xa::Float64,xb::Float64,n::Int)

    xshift = 0.0
    grid = exp.(collect(range(log(xa+xshift),log(xb+xshift),length=n))).-xshift
    grid[1] = xa ; grid[end] = xb

    return grid, xshift
end


"""
    makeknotd(kmin::Float64,kmax::Float64, n::Int, logshift::Float64)
A Function involved in generating the household asset grid. Several methods available.
"""
function makeknotd(kmin::Float64,kmax::Float64, n::Int, logshift::Float64)

    knots = logspaceshift(kmin,kmax,n,logshift)
    return knots
    
end 

function makeknotd(kmin::Float64,kmax::Float64, n::Int)

    knots,logshift = logspaceshift(kmin,kmax,n,1.0,n/4)
    return knots, logshift
end 


"""
    set_parameters()
Function to construct the parameter structure.
"""
function set_parameters( ; β::Float64 = 0.986,    #discount factor
                           γ::Float64 = 2.0,       #Risk aversion
                           ψ::Float64 = 2.0,      #inverse Frisch elasticity
                           ψ1::Float64 = 1.0,     #labor disutility scaling (only in complete markets)
                           B::Float64 = 5.5,      #Supply of Assets: 1.4 times annual GDP = 5.6 times quarterly GDP
                           μ::Float64 = 1.2,      #Markup
                           θ::Float64 = 0.15,     #Probability that Calvo fairy visits
                           Rbar::Float64 = 1.005, #Target quarterly interest rate
                           nz::Int = 3,           #Number of income grid states
                           ρ::Float64 = 0.96566,            #persistence parameter income process
                           σ::Float64 = 0.01695^0.5,          #Std of income process random comp.
                           tax_weights::Array{Float64,1} = [0.0,0.0,1.0],   #tax weights
                           a_min::Float64  = 0.0,            #borrowing constraint
                           a_max::Float64 =75.0,             #maximum asset holdings in grid
                           x_min::Float64 = 0.001,          #parameter to generate the consumption knots
                           nk::Int = 1000,                  #number of grid points used for quadrature grid
                           nb::Int = 200 )                   #number of consumption polynomial knot points

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
    reshape_c(c::Array{Float64},p)

helper function with two methods to replace par2wide and par2long from the MNS code
"""
function reshape_c(c::Array{Float64,1},par::params)
 @unpack nb,nz = par

 c_wide = reshape(c,(nb,nz))

 return c_wide

end


function reshape_c(c::Array{Float64,2},par::params)
 @unpack nb,nz = par

 c_long = reshape(c,(nb*nz,1))

 return c_long[:,1]

end

