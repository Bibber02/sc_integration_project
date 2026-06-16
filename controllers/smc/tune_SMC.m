xstar = [pi;0;0;0];  ustar = 0;     % the equilibrium you want to stabilise
eps = 1e-6;  A = zeros(4);
for i = 1:4
    dx = zeros(4,1); dx(i) = eps;
    A(:,i) = (nonlinearPlant(xstar+dx,ustar) - nonlinearPlant(xstar-dx,ustar))/(2*eps);
end
B = (nonlinearPlant(xstar,ustar+eps) - nonlinearPlant(xstar,ustar-eps))/(2*eps);

Q = diag([100 200 5 10]);  R = 1;     % design knobs
P = icare(A,B,Q,R);
c = (B.'*P).'                       % paste this into the block

% sanity check: three sliding poles in the LHP, one ~0 (the s-direction)
eig( (eye(4) - B*(c.'/(c.'*B)))*A )