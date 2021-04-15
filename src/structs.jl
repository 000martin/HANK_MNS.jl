#File that provides structures used throughout the code

using Parameters, Roots



"""
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

    return grid, xshift
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
A structure that stores the parameter vector of the model
"""

struct param
    β::Array(Float64,1)     #discount factor
    γ::Float64     #Risk aversion
    ψ::Float64     #inverse Frisch elasticity
    B::Float64     #Supply of Asets
    μ::Float64     #Markup
    θ::Float64     #Probability that Calvo fairy visits
end

