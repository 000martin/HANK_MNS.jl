module HANK_MNS

include("structs.jl")

export  set_parameters, params, steady_state, transition_full

include("steady_state.jl")

export EGM_SS, check_steady_state, get_steady_state

include("transition_full.jl")

export get_transition_full, solve_for_transition

#include("complete_markets.jl")

end

"""
#test code
using Main.HANK_MNS
p0 = set_parameters()

b,y,SS = get_steady_state(p0,0.9,[0.95,0.99])

p1 = set_parameters(β = b)

 tr = get_transition_full(20,200,p1,SS)
"""