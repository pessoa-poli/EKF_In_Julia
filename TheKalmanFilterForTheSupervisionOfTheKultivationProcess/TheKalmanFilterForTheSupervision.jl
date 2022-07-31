using Printf, LinearAlgebra
 
# Parameter values used for the simulation model
Kg = 0.1 # gL⁻¹ Monod Constant Glucose
Ke = 0.1 # gL⁻¹ Monod Constant Ethanol
Ygx = 0.17 # gg⁻¹ Conversion factor glucose to biomass
Yge = 0.46 # gg⁻¹ Conversion factor glucose to Ethanol
Yex = 0.6 # gg⁻¹ Conversion factor ethanol to biomass

# Initial Conditions for the extended Kalman Filter
initial_biomass_concentration = 2.4 # gL⁻¹
initial_glucose_concentration = 5.0 # gL⁻¹
initial_ethanol_concentration = 0.1 # gL⁻¹
initial_max_growth_rate_on_glucose = 0.14 # h⁻¹
initial_max_growth_rate_on_ethanol = 0.07 # h⁻¹
initial_estimation_error_covariance_matrix = Diagonal([0.02 for cat in 1:5]) # 5x5 Diagonal matrix, where the value of the first 3 diag elements is 5*0.02g²L², and the last 2 is 0.02h⁻²

#############################################################################################   
# Symbols for symbolic math calculations
G,E,X,P,t = 0,0,0,0,0
Y_gx,Y_ge,Y_ex,mu1,mu2,K_M_G,K_M_E = 0,0,0,0,0

#Variables and Parameters
initX = [2.5; 6; 0.2; 0.15; 0.08] # Initial State (Biomass, Glucose, Ethanol, MaxGrowthRate on Glucose, Max Growth Rate on Ethanol)
initP = Diagonal([0.1,0.02,0.02,0.02,0.02]) # Initial Process Estimation / Initial Estimation Error Covariance Matrix
(initPRows, initPColumns) = size(initP)
init = [initX;reshape(initP,initPRows * initPColumns)] # Combined initial value vector

# For the odesolver
H = [0,0,1,0,0] # Observation Matrix
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

#Model OD
function myODE!(sol, X,G,E,umaxG, umaxE)
    sol[1]=X*( mue1 + mue2); # Biomass
    sol[2]=X*-mue1/Y_gx; # Glucose
    sol[3]=( mue1/Y_gx*Y_ge - mue2/Y_ex); # Ethanol
    sol[4]=0; # mue1
    sol[5]=0; # mue2
end

