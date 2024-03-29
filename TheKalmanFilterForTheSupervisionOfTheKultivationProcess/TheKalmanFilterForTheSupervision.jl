using Printf
println("Starting program.\n")
println("Initializing imported libraries and files.")
using LinearAlgebra, Symbolics, DifferentialEquations
include("./data/BC2_eth_pred.jl")
include("./data/BC2.jl")
include("./data/ME.jl")
include("./data/time.jl")
include("./data/timeE.jl")
include("./data/P.jl")
# using Plots; gr();
using Plots;

println("Defining symbolic variables, initial value and system matrices...")
# Charts file:
charts_file1 = "./charts/BGE.png"
charts_file2 = "./charts/MGR.png"
# Symbols for symbolic math calculations
@variables G E X t

@variables Y_gx Y_ge Y_ex mu1 mu2 K_M_G K_M_E

#Variables / Parameters
initX = [2.5; 6.0; 0.2; 0.15; 0.08] # Initial State (Biomass, Glucose, Ethanol,
# MaxGrowthRate on Glucose, Max Growth Rate on Ethanol) = (X, G, E, mu1, mu2)
initP = Diagonal([0.1,0.02,0.02,0.2,0.02]) # Initial Process Estimation / Initial Estimation Error Covariance Matrix
#(initPRows, initPColumns) = size(initP)
init = [initX;initP[:]] # Combined initial value vector

# For the odesolver
H = [0 0 1 0 0] # Observation Matrix
Q = Diagonal([0.001,0.001,0.001,0.001,0.001]) # Process Noise Covariance Matrix

R = 0.05 # Measurement noise covariance matrix
K1 = 0.1 # Monod Constant Glucose
K2 = 0.1 # Monod Constant Ethanol


#Estimated Parameter values
Ygx=0.15; # Yield glucose -> biomass
Yge=0.34; # Yield glucose -> ethanol
Yex=0.43; # Yield ethanol -> biomass

#Monod terms
mue1 = mu1*G /(G+K_M_G);
mue2 = mu2*E /(E+K_M_E) * (1 - mue1/mu1);

#Model ODE
dS = X * [ ( mue1  + mue2) ; # Biomass
 -mue1/Y_gx  ; # Glucose
 ( mue1/Y_gx*Y_ge - mue2/Y_ex) ; # Ethanol
 0;   # Maximum Specific Growth Rate on Glucose
 0;   # Maximum Specific Growth Rate on Ethanol
 ];

# Jacobian of Model with respect to state variables
F = Symbolics.jacobian(dS,[X,G,E,mu1,mu2]) # jacobian(matrix with each function separated by ';', list with the dependent variables)
# P = initP

dP = F * P + P * F' + Q # Defines the function dP

# Simulation / State prediction and filtering
# Replace all symbolic parameters with their respective numeric values where needed
F  = substitute(F,  Dict(Y_gx=>Ygx, Y_ge=>Yge, Y_ex=>Yex, K_M_G=>K1, K_M_E=>K2));
dS = substitute(dS, Dict(Y_gx=>Ygx, Y_ge=>Yge, Y_ex=>Yex, K_M_G=>K1, K_M_E=>K2));
dP = substitute(dP, Dict(Y_gx=>Ygx, Y_ge=>Yge, Y_ex=>Yex, K_M_G=>K1, K_M_E=>K2));
dS = append!(dS,dP[:])

# Build the ODE function to feed the solver
function tkftcp!(du,u,p,t) # The Kalman Filter For the Cultivation Process's equation
  # du[1] = substitute(dS[1], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
  # du[2] = substitute(dS[2], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
  # du[3] = substitute(dS[3], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
  # du[4] = substitute(dS[4], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
  # du[5] = substitute(dS[5], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
  # du[4] = u[4]
  # du[5] = u[5]
  substitutionDic = merge(
              Dict(X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5]),
              Dict( zip(P[:],u[6:end]) )
      )

  # substitutionDic = Dict(X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5])


  du .= [substitute(someExpression,substitutionDic).val for someExpression=dS]
end

#### PLOTTING WITHOUT UPDATE
# tspan = (0.0,timeE[end])
# prob = ODEProblem(tkftcp!,initX,tspan)
# solutions = []
# sol = solve(prob,save_everystep=true)
# plot(sol)
# savefig("./charts/result.png")
############################

## Y values to store at each step, for posterior plotting of results
Y_Biomass = []
Y_Glucose = []
Y_Ethanol = []
Y_KMG     = []
Y_KME     = []
Y_P = []
currentTime = 0.0

println("Starting Simulation run&record phase.")
for t_index = eachindex(timeE)
  # break
  # global initX, currentTime,P, sol, i # specify variables that need to be valid in the outer scope.
  global init, currentTime

  i = timeE[t_index]

  tspan = (currentTime,i)
  prob = ODEProblem(tkftcp!, init, tspan)
  sol = solve(prob, save_everystep=false)
  MS = ME[t_index]
  newP = reshape(sol(i)[6:end],5,5); # process covariance matrix
  res1 = (newP*H')
  res2 = H*newP*H'
  s = size(res2)
  res3 = (res2 + ones(s[1],s[2])*R)
  K = res1/res3 # Kalman Gain matrix
  PS = sol(i)[1:5]
  FS = PS + K*(MS-PS[3]) #Filtered state
  Pfilt = newP-K*H*newP # filtered process covariance matrix

  append!(Y_Biomass, FS[1])
  append!(Y_Glucose, FS[2])
  append!(Y_Ethanol, FS[3])
  append!(Y_KMG,     FS[4])
  append!(Y_KME,     FS[5])
  append!(Y_P, Pfilt)

  # Update initial state variables for the next iteration.
  currentTime = tspan[2]
  # init = [FS; reshape(newP,25,1)]
  init = [FS; Pfilt[:]]

end

# Plotting
println("Starting Plotting Phase")

println("Plotting Biomass, Glicose and Ethanol.")
plot(timeE, Y_Biomass, label="Biomass")
plot!(timeE, Y_Glucose, label="Glucose")
plot!(timeE, Y_Ethanol, label="Ethanol")
savefig(charts_file1)
println("Saved chart to $charts_file1.")
#
println("Plotting MaxGrowthRates on Glicose and Ethanol")
plot(timeE, Y_KMG, label="KM_G")
plot!(timeE, Y_KME, label="KM_E")
# plot!(timeE,ME, label="Real Ethanol Measurement")
savefig(charts_file2)
println("Saved chart to $charts_file2.")
println("Fin")
