# Complete markets transition path


# Solve for steady state
A = 1
ppi = 1
wage = 1/params.mu
Y = (wage/params.psi1)^(1/(params.γ + params.ψ))
R = 1/params.beta 
C = Y
S = 1
N = S*Y/A
L = N
pbarB = Y/(1- params.beta * (1-params.theta))
pbarA = wage*params.mu*Y/(1-params.beta * (1-params.theta))
pstar = 1
dividend = Y - wage*N

name = ["R","A", "S", "wage","dividend","ppi","Y","C","N","L","pbarB","pbarA","pstar"]

nvar = length(name)
stst = zeros(nvar,1)

for i_ = 1:nvar
    eval(['stst(i_) = ' names{i_} ''])
    eval(['ind_' names{i_} ' = i_'])
end

# Helper functions

function invdistr(Pi)
    @assert abs.(sum.(Pi)-1)<1e-10
    opts.disp = 0
    eigvals

    return
end



function SolveForTransitionCompMkts(Rpath, wagepath, dividendpath, Spath, stst, names, p::params=p)
    
    nvar = length(names)

    for i in eachindex(nvar)
        eval(["ind_' names[i_] ' = i_ "])
    end

    T = length(Rpath)

    X = repeat(stst,1,T)

    for outer_it in 1:100

        for inner_it in 1:100
            
            for t in T-1:-1:2
                #solve backwards using Euler
                X(ind_C,t) = (Params.beta * Rpath(t))^(-1/Params.sigma) * X(ind_C,t+1)
            end
    
            X(ind_Y,2:T-1) = X(ind_C,2:T-1)
            X(ind_N,2:T-1) = Spath(2:T-1).*X(ind_Y,2:T-1)
    
            #labor supply C^(-sigma) * w = psi1 * L^Params.psi2
            X(ind_L,2:T-1) = (wagepath(2:T-1) .* X(ind_C,2:T-1).^(-Params.sigma) / Params.psi1).^(1/Params.psi2)
                
    
            plot(X([ind_N ind_L],:)')
            drawnow
            
            
            #adjust wage, dividend
            oldwage = wagepath
            wagepath(2:T-1) = wagepath(2:T-1).*(X(ind_N,2:T-1)./X(ind_L,2:T-1)).^Params.psi2
            dividendpath(2:T-1) = (X(ind_Y,2:T-1) - wagepath(2:T-1).*X(ind_N,2:T-1))
            
            test = max(abs(wagepath./oldwage - 1))
            disp(['Residual in wage path: ' num2str(test)])
            if test < 1e-6, break end

        end
        #Now compute inflation and S
        #solve backwards from the end for pi
        #then solve forwards from the beginning for S
        
        # initialize these with steady state values
        pbarA = stst(ind_pbarA)
        pbarB = stst(ind_pbarB)
        
        # solve for inflation
        pstar= ones(1,T-2)
        for t = T-1:-1:2
            pbarA = Params.mu*wagepath(t) * X(ind_Y,t) + Params.beta*(1-Params.theta) * X(ind_ppi,t+1)^(- Params.mu/(1- Params.mu)) *pbarA
            pbarB = X(ind_Y,t) + Params.beta*(1-Params.theta)* X(ind_ppi,t+1)^(- 1/(1- Params.mu))*pbarB
            pstar(t) = pbarA/pbarB
            X(ind_ppi,t) = ((1-Params.theta)/(1-Params.theta*pstar(t)^(1/(1-Params.mu))))^(1-Params.mu)
        end
        
        # solve for S
        oldS = Spath
        Spath = ones(1,T)
        Slast = 1
        for t = 2:T-1
            Spath(t) = (1-Params.theta)*Slast*X(ind_ppi,t)^(-Params.mu/(1-Params.mu)) + Params.theta*pstar(t)^(Params.mu/(1-Params.mu))
            Slast = Spath(t)
        end
        
        test = max(abs(Spath./oldS - 1 ))
        disp(['Residual in S path: ' num2str(test)])
        if test < 1e-6, break end
    end
    
    X(ind_wage,:) = wagepath
    eqm = X
    
end
    
