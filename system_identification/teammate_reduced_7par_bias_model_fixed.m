function [dx, y] = teammate_reduced_7par_bias_model_fixed(t, x, u, ...
    p_alpha, p_g1, p_g2, p_b1, p_b2, p_u, p_0, aux, varargin)
%TEAMMATE_REDUCED_7PAR_BIAS_MODEL
% Seven-parameter reduced grey-box model based on the report coordinates.
%
% This is the teammate six-parameter model plus one constant normalized
% torque offset p_0 at the actuated joint.
%
% Estimated parameters:
%   p_alpha = alpha / beta
%   p_g1    = ((m1*c1 + m2*l1)*g) / beta
%   p_g2    = (m2*c2*g) / beta
%   p_b1    = b1_star / beta, with motor back-EMF lumped into b1_star
%   p_b2    = b2 / beta
%   p_u     = k_u / beta
%   p_0     = tau_0 / beta
%
% Not estimated independently:
%   p_c = gamma / beta = (l1/g) * p_g2
%
% States are local measured coordinates if offsets were removed in the
% runner:
%   x(1) = theta_1_local [rad]
%   x(2) = theta_2_local [rad], relative passive joint angle
%   x(3) = omega_1       [rad/s]
%   x(4) = omega_2       [rad/s]
%
% The physical angles used inside sin/cos terms are restored as:
%   theta_1_phys = theta_1_local + aux.theta1_offset
%   theta_2_phys = theta_2_local + aux.theta2_offset
%
% If offsets are not removed in the runner, aux.theta1_offset and
% aux.theta2_offset are zero, so local and physical angles are identical.
%
% Constant torque offset is included through p_0. No actuator lag,
% no Coulomb friction and no Stribeck friction are included in this model.

%#ok<*INUSD>

% In some MATLAB/idnlgrey versions FileArgument arrives as a 1x1 cell
% containing the struct. Unpack it before using dot indexing.
if nargin < 11 || isempty(aux)
    error(['Missing auxiliary constants. Set model0.FileArgument = {aux}, ', ...
           'where aux has fields l1, g, theta1_offset and theta2_offset.']);
end

if iscell(aux)
    if numel(aux) == 1 && isstruct(aux{1})
        aux = aux{1};
    else
        error('Expected aux to be a struct or a 1x1 cell containing a struct.');
    end
end

requiredAuxFields = {'l1', 'g', 'theta1_offset', 'theta2_offset'};
for kAux = 1:numel(requiredAuxFields)
    if ~isfield(aux, requiredAuxFields{kAux})
        error('Missing aux.%s in FileArgument.', requiredAuxFields{kAux});
    end
end

theta1_local = x(1);
theta2_local = x(2);
omega1       = x(3);
omega2       = x(4);

if isempty(u)
    u_in = 0;
else
    u_in = u(1);
end

% Known constants and sensor offsets supplied by the runner.
l1 = aux.l1;
g_const = aux.g;
theta1_offset = aux.theta1_offset;
theta2_offset = aux.theta2_offset;

% Restore physical angles for the nonlinear mechanics.
theta1 = theta1_local + theta1_offset;
theta2 = theta2_local + theta2_offset;

% Teammate's reduction: p_c is fixed by p_g2 and known l1, g.
p_c = (l1 / g_const) * p_g2;

c = cos(theta2);
s = sin(theta2);

% Normalised inertia matrix.
M11 = p_alpha + 1 + 2*p_c*c;
M12 = 1 + p_c*c;
M22 = 1;

% Right-hand side after moving non-acceleration terms to the RHS.
rhs1 = p_u*u_in + p_0 ...
     + p_c*s*(2*omega1*omega2 + omega2^2) ...
     - p_b1*omega1 ...
     + p_g1*sin(theta1) ...
     + p_g2*sin(theta1 + theta2);

rhs2 = -p_c*s*omega1^2 ...
     - p_b2*omega2 ...
     + p_g2*sin(theta1 + theta2);

% Solve M*qddot = rhs.
detM = M11*M22 - M12*M12;

% Numerical guard only. A physically consistent fit should keep detM > 0.
if abs(detM) < 1e-10
    detM = sign(detM + eps) * 1e-10;
end

theta1_ddot = ( M22*rhs1 - M12*rhs2) / detM;
theta2_ddot = (-M12*rhs1 + M11*rhs2) / detM;

dx = zeros(4, 1);
dx(1) = omega1;
dx(2) = omega2;
dx(3) = theta1_ddot;
dx(4) = theta2_ddot;

% Outputs are the same coordinates used in the identification data.
y = zeros(2, 1);
y(1) = theta1_local;
y(2) = theta2_local;

end
