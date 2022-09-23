using Printf, LinearAlgebra, Symbolics, DifferentialEquations
include("./TheKalmanFilterForTheSupervisionOfTheKultivationProcess/data/BC2_eth_pred.jl")
include("./TheKalmanFilterForTheSupervisionOfTheKultivationProcess/data/BC2.jl")
include("./TheKalmanFilterForTheSupervisionOfTheKultivationProcess/data/ME.jl")
include("./TheKalmanFilterForTheSupervisionOfTheKultivationProcess/data/time.jl")
include("./TheKalmanFilterForTheSupervisionOfTheKultivationProcess/data/timeE.jl")
include("./TheKalmanFilterForTheSupervisionOfTheKultivationProcess/data/P.jl")
using Plots

# Symbols for symbolic math calculations
@variables G E X P t
@variables Y_gx Y_ge Y_ex mu1 mu2 K_M_G K_M_E

#Variables and Parameters
initX = [2.5; 6.0; 0.2; 0.15; 0.08] # Initial State (Biomass, Glucose, Ethanol,
  # MaxGrowthRate on Glucose, Max Growth Rate on Ethanol) = (X, G, E, K_M_G, K_M_E)
initP = Diagonal([0.1,0.02,0.02,0.02,0.02]) # Initial Process Estimation / Initial Estimation Error Covariance 
init = [initX; initP[:]]

# For the odesolver
H = [0 0 1 0 0] # Observation Matrix
Q = Diagonal([0.001,0.001,0.001,0.001,0.001]) # Process Noise Covariance Matrix

# R = ones(5,5) * 0.5 # Measurement noise covariance matrix
R = 0.5 # Measurement noise covariance matrix
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

# P is initialized as a matrix of symbolics in a separate file:
# ./TheKalmanFilterForTheSupervisionOfTheKultivationProcess/data/P.jl
dP = F * P + P * F' + Q # Defines the function dP

# Simulation / State prediction and filtering
# Replace all symbolic parameters with their respective numeric values where needed
F = substitute(F, Dict(Y_gx=>Ygx, Y_ge=>Yge, Y_ex=>Yex, K_M_G=>K1, K_M_E=>K2));
dS = substitute(dS, Dict(Y_gx=>Ygx, Y_ge=>Yge, Y_ex=>Yex, K_M_G=>K1, K_M_E=>K2));
dP = substitute(dP, Dict(Y_gx=>Ygx, Y_ge=>Yge, Y_ex=>Yex, K_M_G=>K1, K_M_E=>K2));

println("Printing dP's size")
display(size(dP))

# Building init vector to feed the ODESolver
(initPRows, initPColumns) = size(initP)
init = [initX; reshape(initP,initPRows * initPColumns)] # Combined initial value vector
println("Analyzing init")
println(size(init)[1])
println("fin")
exit(0)

arr=[dS dP]
dSdPSize = size(arr)
println("[dS dP] size = $dSdPSize")
u=init
arr2 = [substitute(arr[cat], Dict([ X=>u[1], G=>u[2], E=>u[3], mu1=>u[4], mu2=>u[5]])) for cat=[1:size(arr)[1]]]
println("Arr2 is: $arr2")
someu = size(arr2)
println("Size of someu is: $someu")