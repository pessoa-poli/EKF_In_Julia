using Printf, LinearAlgebra, Symbolics, DifferentialEquations

include("./data/BC2_eth_pred.jl")
include("./data/BC2.jl")
include("./data/ME.jl")
include("./data/time.jl")
include("./data/timeE.jl")
# include("./data/P.jl")
using Plots

# Charts file:
charts_file = "./charts/result_chart_with_update.png"
# Symbols for symbolic math calculations
@variables G E X P t

@variables Y_gx Y_ge Y_ex mu1 mu2 K_M_G K_M_E

#Variables and Parameters
initX = [2.5; 6.0; 0.2; 0.15; 0.08] # Initial State (Biomass, Glucose, Ethanol,
  # MaxGrowthRate on Glucose, Max Growth Rate on Ethanol) = (X, G, E, K_M_G, K_M_E)
initP = Diagonal([0.1,0.02,0.02,0.02,0.02]) # Initial Process Estimation / Initial Estimation Error Covariance Matrix
(initPRows, initPColumns) = size(initP)
init = [initX;reshape(initP,initPRows * initPColumns)] # Combined initial value vector

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
P = initP
dP = F * P + P * F' + Q # Defines the function dP

# Simulation / State prediction and filtering
# Replace all symbolic parameters with their respective numeric values where needed
F = substitute(F, Dict(Y_gx=>Ygx, Y_ge=>Yge, Y_ex=>Yex, K_M_G=>K1, K_M_E=>K2));
dS = substitute(dS, Dict(Y_gx=>Ygx, Y_ge=>Yge, Y_ex=>Yex, K_M_G=>K1, K_M_E=>K2));
dP = substitute(dP, Dict(Y_gx=>Ygx, Y_ge=>Yge, Y_ex=>Yex, K_M_G=>K1, K_M_E=>K2));

# function P_Update(oldP, )

# println("Printing dS")
# println(dS)

# Build the ODE function to feed the solver
function tkftcp!(du,u,p,t) # The Kalman Filter For the Cultivation Process's equation
  du[1] = substitute(dS[1], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
  du[2] = substitute(dS[2], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
  du[3] = substitute(dS[3], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
  du[4] = substitute(dS[4], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
  du[5] = substitute(dS[5], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5] ]) ).val
end

#### PLOTTING WITHOUT UPDATE 
# tspan = (0.0,time[end])
# prob = ODEProblem(tkftcp!,initX,tspan)
# solutions = []
# sol = solve(prob,save_everystep=true)
# plot(sol)
# savefig("./charts/result.png")
############################

## Y values to store at each step
Y_Biomass = []
Y_Glucose = []
Y_Ethanol = []
currentTime = 0.0

println("Starting Simulation run&record phase.")
for t_index = eachindex(timeE)
  i = timeE[t_index]
  global initX, currentTime,P
  tspan = (currentTime,i)
  prob = ODEProblem(tkftcp!,initX,tspan)
  sol = solve(prob, save_everystep=false)
  currentTime = tspan[2]
  MS = ME[t_index]
  res1 = (P*H')
  res2 = H*P*H'
  s = size(res2)
  res3 = (H*P*H' + ones(s[1],s[2])*R)
  K = res1/res3 # Kalman Gain matrix
  PS = [cat for cat=sol(i)[1:5]]
  FS = PS + K * (MS-PS[3]) #Filtered state
  Pfilt = P-K*H*P # filtered process covariance matrix
  
  # append!(Y_Biomass, sol(i)[1])
  # append!(Y_Glucose, sol(i)[2])
  # append!(Y_Ethanol, sol(i)[3])

  append!(Y_Biomass, FS[1])
  append!(Y_Glucose, FS[2])
  append!(Y_Ethanol, FS[3])

  initX = [sol(i)[1], sol(i)[2], sol(i)[3], sol(i)[4], sol(i)[5]]
  P = Pfilt
end

# Plotting
println("Starting Plotting Phase")
plot(timeE, Y_Biomass, label="Biomass")
plot!(timeE, Y_Glucose, label="Glucose")
plot!(timeE, Y_Ethanol, label="Ethanol")
savefig(charts_file)
println("Saved chart to $charts_file.")
println("Fin")