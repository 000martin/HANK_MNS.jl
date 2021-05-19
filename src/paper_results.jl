
using TypedTables, Plots

"""
    get_figures_3_4(;TR::Int64 = 20, T::Int64 = 200, RChange::Float64 = -0.005)

Function that replicates Figures 3 and 4 from the MNS paper.
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

p1 = plot(0:2*TR-1,tp.Y[1:2*TR]./SS.Y .- 1.0, title="Figure 3: Output response",xlabel = "Quarter", ylabel = "Output")  
#plot!([TR+1],seriestype = "vline", label = "Time of Rate Change")
p2 = plot(0:2*TR-1,tp.pΠ[1:2*TR] .- 1.0,title = "Figure 4: Inflation Response",xlabel = "Quarter", ylabel = "Inflation")
#plot!([TR+1],seriestype = "vline", label = "Time of Rate Change")
plot(p1,p2, layout = (2,1), legend = false)
end

"""
    get_figures_5_6(; Horizon::StepRange{Int64,Int64} = 1:2:41, T::Int64 = 200, RChange::Float64 = -0.005)

Function to replicate Figures 5 and 6 from the MNS paper.
WARNING: Running this take quite long, as a large number of transition paths are computed
"""
function get_figures_5_6(; Horizon::StepRange{Int64,Int64} = 1:2:41, T::Int64 = 200, RChange::Float64 = -0.005)

    TRs = collect(Horizon)

    #pre-allocate vectors for results
    Y_response = Array{Float64,1}(undef,length(TRs)); Π_response = Array{Float64,1}(undef,length(TRs))

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
    println("Initial Output response: ",tr.Y[1]/SS.Y - 1.0,"Initial Inflation response: ",tr.pΠ[1]-1.0)
    println("")

    end

    p1 = plot(Horizon, Y_response/SS.Y .- 1.0, title="Figure 5: Initial Output response",xlabel = "Horizon rate change", ylabel = "Output")  
    #plot!([TR+1],seriestype = "vline", label = "Time of Rate Change")
    p2 = plot(Horizon,Π_response .- 1.0,title = "Figure 6: Initial Inflation response",xlabel = "Horizon rate change", ylabel = "Inflation")
    #plot!([TR+1],seriestype = "vline", label = "Time of Rate Change")
    plot(p1,p2, layout = (2,1), legend = false)
end

"""
    get_table_2(;Horizon::Int=20,T::Int=200,RChange::Float64 = -0.005)

Function that generates the results in table 2 from the MNS paper.
Returns a TypedTable.jl table containing the results.
Baseline case: Interest rate change 20 periods ahead (Horizon), complete transition period length T=200, 
one time interest rate change (RChange) of 50 basis points
"""

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

 return Table(Case = cases,  initial_output_response = resp_output, initial_inflation_response= resp_inflation)

end


