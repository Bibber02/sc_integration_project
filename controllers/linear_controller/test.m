% 1. Load the saved grey-box model objects from your ID script
load('system_identification/active_rod_identification_result.mat'); 

% The variable 'model_est' is now in your workspace. 
% This is your fully trained idnlgrey object.

% 2. Set the operating point to the upright equilibrium position
% (theta_1 = 0, theta_2 = 0, omega_1 = 0, omega_2 = 0)
model_est.InitialStates(1).Value = pi; 
model_est.InitialStates(2).Value = 0; 
model_est.InitialStates(3).Value = 0; 
model_est.InitialStates(4).Value = 0; 

% 3. Extract the continuous-time linear State-Space model
% The ss() command automatically computes the Jacobians of your ODEs
% using your trained parameters around the initial states we just set.
sys_linear = linearize(model_est);

% 4. Pull out the matrices!
[A, B, C, D] = ssdata(sys_linear);

% Display them in the command window to verify
disp('A Matrix:');
disp(A);
disp('B Matrix:');
disp(B);
