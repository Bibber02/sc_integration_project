function [u, e1dot] = smc(x, r, dr, z)      %#codegen
% x  : 4x1 state estimate [th1; th2; dth1; dth2] from EKF
% r  : theta1 reference                (scalar)
% dr : theta1 reference rate           (scalar; 0 for step refs)
% z  : integral state from Simulink integrator  (= integral of e1)
    x = x(:);

% ---- parameters (must match identified plant) ----
    pa=64.7915; pb1=3.2117e3; pc1=635.3865; pg1=430.5226; pu=2.2424e4;
    pb2=0.0670; pg2=112.0268; pc2=0.2491; psdelta2=0.1836;
    epsv1=0.05; epsv2=0.03; vs2=3.5;  l1=0.1; g=9.81;
    pc=(l1/g)*pg2;

% ---- design constants (from augmented icare) ----
    c   = [ ... ; ... ; ... ; ... ];   % paste ca(1:4)
    ci  = ...;                          % paste ca(5)
    k   = 30;        % switching gain  -- TUNE
    phi = 0.3;       % boundary layer  -- TUNE (integrator now handles SS error)

% ---- tracking error: theta1 -> r, upper link inverted ----
    e    = x - [r; pi; dr; 0];
    e(1) = mod(e(1)+pi,2*pi)-pi;
    e(2) = mod(e(2)+pi,2*pi)-pi;

% ---- augmented sliding variable ----
    s    = c.'*e + ci*z;

% ---- inlined drift f0 = f(x,0) ----
    th1=x(1); th2=x(2); d1=x(3); d2=x(4);
    c2=cos(th2); s2=sin(th2);
    M11=pa+1+2*pc*c2; M12=1+pc*c2; M22=1; detM=M11*M22-M12*M12;
    Cq1=-pc*s2*(2*d1*d2+d2^2);
    Cq2= pc*s2*d1^2;
    F1 = pb1*d1 + pc1*tanh(d1/epsv1);
    strib = pc2 + psdelta2*exp(-(d2/vs2)^2);
    F2 = pb2*d2 + strib*tanh(d2/epsv2);
    G1 = -pg1*sin(th1) - pg2*sin(th1+th2);
    G2 = -pg2*sin(th1+th2);
    r1 = -Cq1 - F1 - G1;
    r2 = -Cq2 - F2 - G2;
    dd1 = ( M22*r1 - M12*r2)/detM;
    dd2 = (-M12*r1 + M11*r2)/detM;
    f0  = [d1; d2; dd1; dd2];

% ---- input gain c'*g(x), guarded against the sign change ----
    cg  = (pu/detM)*( c(3) - c(4)*M12 );
    cgmin = 1e-3;
    if abs(cg) < cgmin
        if cg >= 0, cg = cgmin; else, cg = -cgmin; end
    end

% ---- equivalent control: cancel drift + reference motion + integral term ----
%   sdot = c'(f0 + g u) - c'*xrefdot + ci*e1 = 0,  xrefdot ~ [dr;0;0;0]
    ueq = -( c.'*f0 - c(1)*dr + ci*e(1) )/cg;

% ---- switching term ----
    usw = -(k/cg)*max(min(s/phi,1),-1);

    uun = ueq + usw;
    u   = max(min(uun, 1), -1);

% ---- conditional integration (anti-windup): freeze integrator when saturated ----
    e1dot = e(1);
    if (uun > 1) || (uun < -1)
        e1dot = 0;
    end
end