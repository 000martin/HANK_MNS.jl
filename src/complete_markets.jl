using LinearAlgebra

"""
    get_transition_CompMkts(TR::Int, T::Int, p::params; RChange::Float64 = -0.005)

Solves for the complete markets transition path.
Inputs: TR (time until single-quarter interest rate change), T (time horizon transition path), 
(Change of R at time TR))
"""
function get_transition_CompMkts(TR::Int, T::Int, p::params; RChange::Float64 = -0.005)
    
    # Update struct for representative agent case
    p.ψ1 = p.ψ1 * dot(p.Γ, p.z.^(1+1/p.ψ))
    p.β = 1/p.Rbar
    
    # Solve for steady state
    A = 1
    ppi = 1
    wage = 1/p.μ
    Y = (wage/p.ψ1)^(1/(p.γ + p.ψ))
    R = 1/p.β
    C = Y
    S = 1
    N = S*Y/A
    L = N
    pbarB = Y/(1-p.β * (1-p.θ))
    pbarA = wage*p.μ*Y/(1-p.β * (1-p.θ))
    pstar = 1
    dividend = Y - wage*N

    name = ["R","A", "S", "wage","dividend","ppi","Y","C","N","L","pbarB","pbarA","pstar"]

    nvar = length(name)
    stst = zeros(nvar,1)

    stst[1] = R; stst[2] = A; stst[3] = S; stst[4] = wage; stst[5] = dividend; stst[6] = ppi; stst[7] = Y;
    stst[8] = C; stst[9] = N; stst[10] = L; stst[11] = pbarB; stst[12] = pbarA; stst[13] = pstar 

    ind_R = 1; ind_A = 2; ind_S = 3; ind_wage = 4; ind_dividend = 5; ind_ppi = 6; ind_Y = 7; ind_C = 8; 
    ind_N = 9; ind_L = 10; ind_pbarB = 11; ind_pbarA = 12; ind_pstar = 13
    
    # Initial guesses
    w = 1/p.μ
    R = p.Rbar

    wagepath = fill(w,1,T)
    dividendpath = fill(stst[ind_dividend],1,T)
    Spath = ones(1,T)
    
    # Policy
    Rpath =  fill(R,1,T)
    Rpath[TR+2] =  p.Rbar + RChange
    
    tp = solve_for_transition_CompMkts(Rpath, wagepath, dividendpath, Spath, stst, name, p::params)

    return tp
end


"""
    solve_for_transition_CompMkts(Rpath, wagepath, dividendpath, Spath, stst, name, p::params)

Solves for the complete markets transition path for a given interest rate path, wpage path, dividend path, etc. 
"""
function solve_for_transition_CompMkts(Rpath, wagepath, dividendpath, Spath, stst, name, p::params)
    
    nvar = length(name)

    ind_R = 1; ind_A = 2; ind_S = 3; ind_wage = 4; ind_dividend = 5; ind_ppi = 6; ind_Y = 7; ind_C = 8; 
    ind_N = 9; ind_L = 10; ind_pbarB = 11; ind_pbarA = 12; ind_pstar = 13

    T = length(Rpath)

    X = repeat(stst,1,T)

    for outer_it = 1:100

        for inner_it = 1:100
            
            for t = T-1:-1:2
                # Solve backwards using Euler
                X[ind_C,t] = (p.β * Rpath[t])^(-1/p.γ) * X[ind_C,t+1]
            end
    
            X[ind_Y,2:T-1] = X[ind_C,2:T-1]
            X[ind_N,2:T-1] = Spath[2:T-1].*X[ind_Y,2:T-1]
    
            X[ind_L,2:T-1] = (wagepath[2:T-1] .* X[ind_C,2:T-1].^(-p.γ) / p.ψ1).^(1/p.ψ)
                
            # Adjust wage, dividend
            oldwage = wagepath[:]
            wagepath[2:T-1] = wagepath[2:T-1].*(X[ind_N,2:T-1]./X[ind_L,2:T-1]).^p.ψ
            dividendpath[2:T-1] = (X[ind_Y,2:T-1] - wagepath[2:T-1].*X[ind_N,2:T-1])
         
            test = maximum(abs.(wagepath[2:T-1]./oldwage[2:T-1] .- 1.0))
            println(string("Residual in wage path: ", test))
            if test < 1e-6 break end
        end
        
        # Now compute inflation and S
        # Solve backwards from the end for pi
        # Then solve forwards from the beginning for S
        
        # Initialize these with steady state values
        pbarA = stst[ind_pbarA]
        pbarB = stst[ind_pbarB]
        
        # Solve for inflation
        pstar = ones(1,T-1)
        for t = T-1:-1:2
            pbarA = p.μ * wagepath[t] * X[ind_Y,t] + p.β * (1-p.θ) * X[ind_ppi,t+1]^(-p.μ /(1-p.μ)) * pbarA
            pbarB = X[ind_Y,t] + p.β * (1-p.θ) * X[ind_ppi,t+1]^ (-1/(1-p.μ))*pbarB
            pstar[t] = pbarA/pbarB
            X[ind_ppi,t] = ((1-p.θ)/(1-p.θ*pstar[t]^(1/(1-p.μ))))^(1-p.μ)
        end
        
        # Solve for S
        oldS = Spath
        Spath = ones(1,T)
        Slast = 1

        for t = 2:T-1
            Spath[t] = (1-p.θ)*Slast*X[ind_ppi,t]^(-p.μ/(1-p.μ)) + p.θ*pstar[t]^(p.μ/(1-p.μ))
            Slast = Spath[t]
        end
    
        test = findmax(abs.(Spath./oldS .- 1 ))[1]
        println(string("Residual in S path: ", test))
        if test < 1e-6 break end
    end
    
    X[ind_wage,:] = wagepath
    eqm = X

    Spath = vec(Spath)
    wagepath = vec(wagepath)
    Rpath = vec(Rpath)
    dividendpath = vec(dividendpath)
    Ypath = eqm[ind_Y,2:T-1]./stst[ind_Y].-1
    pΠpath = eqm[ind_ppi,2:T-1].-1

    return transition_CompMkts(Spath,wagepath,pΠpath,Ypath,Rpath,dividendpath)
end