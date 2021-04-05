+module HANK_MNS
export hello, domath, domath2

"""
    hello(who::String)

Return "Hello, `who`".
"""
hello(who::String) = "Hello, $who"

"""
    domath(x::Number)

Return `x + 5`.
"""
domath(x::Number) = x + 5


"""
    domath2(x::Number)
Return `x+6`. Serves as a text.
"""
domath2(x::Number) = x+6
end
