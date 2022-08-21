using Symbolics
using Latexify
using DifferentialEquations
using Plots
gr()

@variables s s0 v0 a t

function s!(du,u,p,t) 
    f = s0 + v0*t + (a*t^2)/2
    du = substitute(f, Dict([a=>du[1], v0=>du[2], s0=>du[3]]))
end


tspan = [1.0, 6.0]
v0 = 10.0
a = 9.8
s0 = 0.0
du = [a,v0,s0]

prob = ODEProblem(s!,du,tspan)
sol = solve(prob)
plot(sol)
savefig("./charts/test4Plot.png")
