
using TypedTables, Plots

"""
    get_figures_3_4(;TR::Int64 = 20, T::Int64 = 200, RChange::Float64 = -0.005)

Convenience function that replicates Figures 3 and 4 from the MNS paper.

Just running `get_figures_3_4()` will display equivalents to Figures 3 and 4 in the MNS paper.

The function relies on `set_parameters`, `get_steady_state`, `get_transition_full` and `get_transition_CompMkts`.

Optional argument: Time horizon of interest rate change (TR). Size of interest rate change (RChange)
"""
function get_figures_3_4(;TR::Int64 = 20, T::Int64 = 200, RChange::Float64 = -0.005)

p = set_parameters()

println("")
println("Solving for Steady State")
println("")
b,y,SS = get_steady_state(p,0.6,[0.95,0.99])

p.β = b

println("")
println("Solving for transition path")
println("")
tp = get_transition_full(TR,T,p,SS; RChange=RChange)

println("")
println("Solving for complete markets transition path")
println("")

tpc = get_transition_CompMkts(TR,T,p; RChange=RChange)

p1 = plot(0:2*TR-1,(tp.Y[1:2*TR]./SS.Y .- 1.0)*100.0, title="Figure 3: Output response",xlabel = "Quarter", ylabel = "Output",
        label = "Incomplete Markets")  
plot!(0:2*TR-1,tpc.Y[1:2*TR]*100.0,label = ["Complete Markets"] )

p2 = plot(0:2*TR-1,(tp.pΠ[1:2*TR] .- 1.0)*100.0,title = "Figure 4: Inflation Response",xlabel = "Quarter", ylabel = "Inflation",
        label = "Incomplete Markets")
plot!(0:2*TR-1,tpc.pΠ[1:2*TR]*100.0,label = ["Complete Markets"] )

return plot(p1,p2, layout = (2,1))
end

"""
    get_figures_5_6(; Horizon::StepRange{Int64,Int64} = 1:2:41, T::Int64 = 200, RChange::Float64 = -0.005)

Convenience function that replicates Figures 5 and 6 from the MNS paper.

Just running `get_figures_5_6()` will display equivalents to Figures 5 and 6 in the MNS paper.

The function relies on `set_parameters`, `get_steady_state`, `get_transition_full` and `get_transition_CompMkts`.
        
Optional argument: 
..* Horizons of interest rate changes to consider (`Horizon`), to be supplied as `StepRange`.
..* Total length of transition period to compute (`T`), i.e. assuming that after T periods economy will be back in SS
..* Size of interest rate change (RChange).

All default values corresppond to the values used by MNS.
"""
function get_figures_5_6(; Horizon::StepRange{Int64,Int64} = 1:2:41, T::Int64 = 200, RChange::Float64 = -0.005)

    TRs = collect(Horizon)

    #pre-allocate vectors for results
    Y_response = Array{Float64,1}(undef,length(TRs)); Π_response = Array{Float64,1}(undef,length(TRs))
    Y_response_CM = Array{Float64,1}(undef,length(TRs)); Π_response_CM = Array{Float64,1}(undef,length(TRs))

    #get initial steady state
    p = set_parameters()
    b, y , SS = get_steady_state(p,0.6,[0.95,0.99])
    p.β = b #update beta


    for horiz = 1:length(TRs)


    println("")
    println("Solving for transition path with horizon ",TRs[horiz])
    println("")

    tr = get_transition_full(TRs[horiz],T,p,SS; RChange = RChange)

    Y_response[horiz] = tr.Y[1]
    Π_response[horiz] = tr.pΠ[1]
    
    println("")
    println("Initial Output response: ",tr.Y[1]/SS.Y - 1.0," Initial Inflation response: ",tr.pΠ[1]-1.0)
    println("")

    trc = get_transition_CompMkts(TRs[horiz],T,p; RChange=RChange)

    Y_response_CM[horiz] = trc.Y[1]
    Π_response_CM[horiz] = trc.pΠ[1]

    p = set_parameters(β = b,ψ1 = 1.0) #ensure par structure is stable

    end

    p1 = plot(Horizon, (Y_response/SS.Y .- 1.0)*100.0, title="Figure 5: Initial Output response",xlabel = "Horizon rate change", ylabel = "Output",
                label = "Incomplete Markets")  
    plot!(Horizon, Y_response_CM*100.0, label = "Complete Markets")
    
    p2 = plot(Horizon,(Π_response .- 1.0)*100.0,title = "Figure 6: Initial Inflation response",xlabel = "Horizon rate change", ylabel = "Inflation",
                label = "Complete Markets")
    plot!(Horizon, Π_response_CM*100.0, label = "Complete Markets")

    return plot(p1,p2, layout = (2,1))
end

"""
    get_table_2(;Horizon::Int=20,T::Int=200,RChange::Float64 = -0.005)

Convenience function that replicates Table 2 from the MNS paper.

Just running `get_table_2()` will return an equivalent to Table 2 as a table object.

The function relies on `set_parameters`, `get_steady_state`, `get_transition_full` and `get_transition_CompMkts`.
        
Optional argument: 
..* Time horizon of interest rate change (`TR`)
..* Total length of transition period to compute (`T`), i.e. assuming that after T periods economy will be back in SS 
..* Size of interest rate change (`RChange`)

Default values correspond to the values chosen by MNS.

### Note: The returned values will not be identical to the ones obtained by MNS. 
### However, notice that chosen unit is Basis Points (=0.01 percent), so the actual numerical difference between the MNS solution and ours is small.
"""
function get_table_2(;Horizon::Int=20,T::Int=200,RChange::Float64 = -0.005)

 #arrays for results
 resp_inflation = Array{Float64,1}(undef,5)
 resp_output    = Array{Float64,1}(undef,5)
 cases = ["Baseline","High Risk","High Asset","High Risk and High Asset","Complete Markets"]

 #output standard deviations for different cases:
 σ2s = [0.01695,0.033,0.01695,0.024]
 #aggregate assets for different cases
 Bs = [5.5,5.5,15.15,15.15]

 #loop over cases
 for case = 1:4

    #adjust guess for β range and range of asset grid depending on case
    if Bs[case] == 15.15
    βmin = 0.97 ;  βmax = 0.995 ; a_max = 110.0 
    else
    βmin = 0.95 ;  βmax = 0.99 ; a_max = 75.0
    end

    #generate initial parameter vector
    p = set_parameters(B = Bs[case], σ = σ2s[case]^0.5, a_max = a_max)
    
    println("")
    println("Solving for Steady State: ",cases[case]," Calibration")
    println("")

    #calibrate β and get steady state
    β,y,SS = get_steady_state(p,0.6,[βmin,βmax])

    p.β = β #update parameter structure

    println("")
    println("Solving for Transition path: ",cases[case]," Calibration")
    println("")

    tr = get_transition_full(Horizon,T,p,SS;RChange = RChange)

    #save results - convert to basis points
    resp_inflation[case] = (tr.pΠ[1] - 1.0)*10000.0 
    resp_output[case]    = (tr.Y[1]/SS.Y - 1.0)*10000.0

 end
 
 p = set_parameters()
 tpc = get_transition_CompMkts(Horizon,T,p;RChange = RChange)
 resp_inflation[5] = tpc.pΠ[1]*10000.0 
 resp_output[5] = tpc.Y[1]*10000.0

 return Table(Case = cases,  initial_output_response = resp_output, initial_inflation_response= resp_inflation)

end


