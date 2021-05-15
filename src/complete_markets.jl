"""
   complete_markets()

"""

function complete_markets()
    
    p = set_parameters()

    # Add twofold definition of ψ to the p struct
    p.ψ1 = 1
    p.ψ2 = 2

    # Update field values for representative agent case
    p.ψ1 = p.ψ1 * dot(p.Γ, p.z.^(1+1/p.ψ2)) 
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

    names = ["R","A", "S", "wage","dividend","ppi","Y","C","N","L","pbarB","pbarA","pstar"]

    nvar = length(names)
    stst = zeros(nvar,1)

    for i in eachindex(nvar)
        eval("stst[i] =" names[i])
        eval("ind_" names[i] " = i")
    end

    # Solve for transition
    T = 200

    # Initial guesses
    [R, w] = steadystateprices()
    wagepath = repeat(w,1,T)
    dividendpath = repeat(stst[ind_dividend],1,T)
    Spath = ones(1,T)
    
    # Policy
    Rpath =  repeat(R,1,T)
    HORIZ = 20
    Rpath[HORIZ+2] = 1
    
    eqm = transition_complete_markets(Rpath, wagepath, dividendpath, Spath, stst, names)
    
end


"""
   transition_complete_markets()

"""

function transition_complete_markets(Rpath, wagepath, dividendpath, Spath, stst, names, p::p)
    
    nvar = length(names)

    for i in eachindex(nvar)
        eval("ind_" names[i] = i)
    end

    T = length(Rpath)

    X = repeat(stst,1,T)

    for outer_it in 1:100

        for inner_it in 1:100
            
            for t = T-1:-1:2
                #solve backwards using Euler
                X[ind_C,t] = (p.β * Rpath[t])^(-1/p.σ) * X[ind_C,t+1]
            end
    
            X[ind_Y,2:T-1] = X[ind_C,2:T-1]
            X[ind_N,2:T-1] = Spath[2:T-1].*X[ind_Y,2:T-1]
    
            #labor supply C^(-sigma) * w = psi1 * L^p.psi2
            X[ind_L,2:T-1] = (wagepath[2:T-1] .* X[ind_C,2:T-1].^(-p.σ) / p.ψ1).^(1/p.ψ2)
                
    
            #plot(X[[ind_N ind_L],:]')
            #drawnow
            
            
            #adjust wage, dividend
            oldwage = wagepath
            wagepath[2:T-1] = wagepath[2:T-1].*(X[ind_N,2:T-1]./X[ind_L,2:T-1]).^p.ψ2
            dividendpath[2:T-1] = (X[ind_Y,2:T-1] - wagepath[2:T-1].*X[ind_N,2:T-1])
            
            test = max(abs(wagepath./oldwage - 1))
            println(string("Residual in wage path: ", test))
            if test < 1e-6, break end
        end
        
        # Now compute inflation and S
        # Solve backwards from the end for pi
        # Then solve forwards from the beginning for S
        
        # Initialize these with steady state values
        pbarA = stst[ind_pbarA]
        pbarB = stst[ind_pbarB]
        
        # Solve for inflation
        pstar = ones(1,T-2)
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
            Spath[t] = (1-p.θ)*Slast*X(ind_ppi,t)^(-p.μ/(1-p.μ)) + p.θ*pstar(t)^(p.μ/(1-p.μ))
            Slast = Spath[t]
        end
        
        test = max(abs(Spath./oldS - 1 ))
        println(string("Residual in S path: ", test))
        if test < 1e-6, break end
    end
    
    X[ind_wage,:] = wagepath
    eqm = X
    
end