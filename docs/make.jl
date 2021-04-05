push!(LOAD_PATH,"../src/")
using Documenter, HANK_MNS

makedocs(modules = [HANK_MNS], sitename = "HANK_MNS.jl")

deploydocs(repo = "github.com/000martin/HANK_MNS.jl.git", devbranch = "main")
