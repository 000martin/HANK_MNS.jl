"""
    get_table_2()

Function that generates the results in table 2 from the MNS paper.
Returns a TypedTable.jl table containing the results.
Baseline case: Interest rate change 20 periods ahead (Horizon), complete transition period length T=200, 
one time interest rate change (RChange) of 50 basis points
"""

using TypedTables

function get_table_2(;Horizon::Int=20,T::Int=200,RChange::Float64 = -0.005)

 #arrays for results
 resp_inflation = Array{Float64,1}(undef,4)
 resp_output    = Array{Float64,1}(undef,4)
 cases = ["Baseline","High Risk","High Asset","High Risk and High Asset"]

 #output standard deviations for different cases:
 σ2s = [0.01695,0.033,0.01695,0.024]
 #aggregate assets for different cases
 Bs = [5.5,5.5,15.15,15.15]

 #loop over cases
 for case = 1:4

    #adjust guess for β range and range of asset grid depending on case
    if B[case] == 15.15
    βmin = 0.97 ;  βmax = 0.995 ; a_max = 110.0 
    else
    βmin = 0.95 ;  βmax = 0.99 ; a_max = 75.0
    end

    #generate initial parameter vector
    p = set_parameters(B = Bs[case], σ = σ2s[case]^0.5, a_max = a_max)

    #calibrate β and get steady state
    β,y,SS = get_steady_state(p,0.6,[βmin,βmax])

    p.β = β #update parameter structure

    tr = get_transition_full(Horizon,T,p,SS;RChange = RChange)

    #save results - convert to basis points
    resp_inflation[case] = (tr.pΠ[1] - 1.0)*10000.0 
    resp_output[case]    = (tr.Y[1]/SS.Y - 1.0)*10000.0

 end

 return Table(Case = cases,  initial_output_response = resp_output, initial_inflation_response= resp_inflation)

end