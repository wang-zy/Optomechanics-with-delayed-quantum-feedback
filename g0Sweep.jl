push!(LOAD_PATH, pwd())
ID = addprocs(16)
using TwoPhotonSweep

g0 = 20.1:0.1:80
pmap(g0Sweep, g0)
rmprocs(ID)
