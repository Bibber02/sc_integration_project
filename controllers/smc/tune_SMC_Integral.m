xstar = [pi; pi; 0; 0];   ustar = 0;     % inverted equilibrium (must match controller)

eps = 1e-6;  A = zeros(4);
for i = 1:4
    dx = zeros(4,1); dx(i) = eps;
    A(:,i) = (nonlinearPlant(xstar+dx,ustar) - nonlinearPlant(xstar-dx,ustar))/(2*eps);
end
B = (nonlinearPlant(xstar,ustar+eps) - nonlinearPlant(xstar,ustar-eps))/(2*eps);

% augment:  zdot = e1  (integral of the theta1 tracking error)
Aa = [A,         zeros(4,1);
      1 0 0 0,   0         ];
Ba = [B; 0];

Qa = diag([100 200 5 10  50]);   R = 1;   % last weight (50) drives the integrator -- TUNE
Pa = icare(Aa, Ba, Qa, R);
ca = (Ba.'*Pa).';                          % 5x1

c  = ca(1:4)                               % paste into the block
ci = ca(5)                                 % paste into the block

% sanity check: four sliding poles in the LHP, one ~0 along s
eig( (eye(5) - Ba*(ca.'/(ca.'*Ba)))*Aa )