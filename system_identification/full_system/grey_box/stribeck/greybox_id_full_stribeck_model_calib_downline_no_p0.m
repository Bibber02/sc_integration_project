function [dx, y] = greybox_id_full_stribeck_model_calib_downline_no_p0(t, x, u, p_a, p_b1, p_c1, p_g1, p_u, p_b2, p_g2, p_c2, p_sdelta2, v_s2, eps_v1, eps_v2, theta_scale, theta1_offset, theta_abs_down_raw, varargin)
%GREYBOX_ID_FULL_STRIBECK_MODEL_CALIB_DOWNLINE_NO_P0
% Reduced full-system grey-box model with sensor calibration correction,
% no constant motor torque bias, and a down-line calibration constraint.
%
% State convention:
%   x(1) = theta_1 in the saved/measured coordinate [rad]
%   x(2) = theta_2 in the saved/measured coordinate [rad]
%   x(3) = theta_1_dot in the saved/measured coordinate [rad/s]
%   x(4) = theta_2_dot in the saved/measured coordinate [rad/s]
%
% Calibration correction on already-calibrated measured angles:
%   theta_1_phys = theta_scale * theta_1_meas + theta1_offset
%   theta_2_phys = theta_scale * theta_2_meas + theta2_offset
%
% The same scale is used for both angles. The theta_2 offset is computed
% from the initial-rest hanging-down line:
%   if theta_1_meas + theta_2_meas = theta_abs_down_raw,
%   then theta_1_phys + theta_2_phys = pi.
%
% Therefore:
%   theta2_offset = pi - theta_scale*theta_abs_down_raw - theta1_offset
%
% theta_abs_down_raw is fixed in the identification script from the mean
% measured theta_1 + theta_2 during the first rest segment before input starts.

%#ok<NASGU> t is required by idnlgrey but not used explicitly.

% Known constants
g  = 9.81;   % [m/s^2]
l1 = 0.10;   % [m]

% Keep calibration scale away from zero for numerical safety.
theta_scale = max(theta_scale, 1e-6);

% Dependent reduced coupling parameter
p_c = l1 * p_g2 / g;

% Computed theta_2 offset from the down-line constraint.
theta2_offset = pi - theta_scale * theta_abs_down_raw - theta1_offset;

% Measured-coordinate states
th1_m  = x(1);
th2_m  = x(2);
dth1_m = x(3);
dth2_m = x(4);

% Corrected physical coordinates used inside the dynamics
th1  = theta_scale * th1_m  + theta1_offset;
th2  = theta_scale * th2_m  + theta2_offset;
dth1 = theta_scale * dth1_m;
dth2 = theta_scale * dth2_m;

% Input voltage
u_1 = u(1);

% Trigonometric terms
cos_th2 = cos(th2);
sin_th2 = sin(th2);

% Joint 1 friction: Coulomb + viscous, evaluated in physical velocity
eps_v1 = max(eps_v1, 1e-5);
F_1 = p_b1 * dth1 + p_c1 * tanh(dth1 / eps_v1);

% Joint 2 friction: Stribeck + Coulomb + viscous, evaluated in physical velocity
v_s2   = max(v_s2, 1e-5);
eps_v2 = max(eps_v2, 1e-5);

friction_level_2 = p_c2 + p_sdelta2 * exp(-(dth2 / v_s2)^2);
F_2 = p_b2 * dth2 + friction_level_2 * tanh(dth2 / eps_v2);

% Reduced inertia matrix M_r(q)
M_11 = p_a + 1 + 2 * p_c * cos_th2;
M_12 = 1 + p_c * cos_th2;
M_21 = M_12;
M_22 = 1;

M = [M_11, M_12;
     M_21, M_22];

% Reduced Coriolis/centrifugal vector C_r(q,qdot) qdot
C_11 = -p_c * sin_th2 * dth2;
C_12 = -p_c * sin_th2 * (dth1 + dth2);
C_21 =  p_c * sin_th2 * dth1;
C_22 = 0;

C_vec = [C_11 * dth1 + C_12 * dth2;
         C_21 * dth1 + C_22 * dth2];

% Reduced gravity vector, evaluated at corrected physical angles
G = [-p_g1 * sin(th1) - p_g2 * sin(th1 + th2);
     -p_g2 * sin(th1 + th2)];

% Input vector: motor voltage gain only, no constant torque bias
I = [p_u * u_1;
     0];

% Dynamics in corrected physical coordinates:
%   M*qddot + C*qdot + F + G = I
rhs = I - C_vec - [F_1; F_2] - G;
qdd = M \ rhs;

% Convert physical accelerations back to the saved/measured coordinate.
dx = zeros(4, 1);
dx(1) = dth1_m;
dx(2) = dth2_m;
dx(3) = qdd(1) / theta_scale;
dx(4) = qdd(2) / theta_scale;

% Output vector: saved/measured coordinate
y = zeros(2, 1);
y(1) = th1_m;
y(2) = th2_m;

end
