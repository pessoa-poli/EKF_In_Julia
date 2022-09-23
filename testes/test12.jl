using DifferentialEquations, Plots, Symbolics
@variables σ,y,x,ρ,z,β
lorenz_matrix = [σ*(y-x);
                 x*(ρ-z)-y;
                 x*y-β*z]

# Build array of callables from the matrix
array_of_symbolic_function_variables = 
callables = [eval(build_function(cat,Symbolics.get_variables(cat))) for cat=lorenz_matrix[:]]
println("Displaying callables.")
display(callables)
exit(0)

function lorenz!(du,u,p,t)
    du = []
    # du[1] = 10.0*(u[2]-u[1])
    # du[2] = u[1]*(28.0-u[3]) - u[2]
    # du[3] = u[1]*u[2] - (8/3)*u[3]
   end

u0 = [1.0,0.0,0.0]
tspan = (0.0,1.0)
p = [10.0,28.0,8/3]
prob = ODEProblem(parameterized_lorenz!,u0,tspan,p)

sol = solve(prob)

plot(sol,vars=(1,2,3))
plot(sol,vars=(0,2))