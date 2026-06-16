function u = smc(x)            %#codegen
    x = x(:);                      % force 4x1 column

% ---- integrator state: integral of the theta1 tracking error ----   NEW
    persistent xi
    if isempty(xi); xi = 0; end
    h = 0.01;                      % sampling period [s] -- MUST equal workspace h   NEW

% ---- parameters ----
    pa=64.7915; pb1=3.2117e3; pc1=635.3865; pg1=430.5226; pu=2.2424e4;
    pb2=0.0670; pg2=112.0268; pc2=0.2491; psdelta2=0.1836;
    epsv1=0.05; epsv2=0.03; vs2=3.5;  l1=0.1; g=9.81;
    pc=(l1/g)*pg2;

% ---- design constants ----
    c      = [19.9432; 27.7084; 0.3566; 1.5661];   % B'P sliding vector
    k      = 30;                   % switching gain   -- TUNE
    phi    = 1.0;                  % boundary layer   -- TUNE
    ki     = 0.5;                  % integral gain    -- TUNE
    xi_max = 2.0;                  % anti-windup clamp
    xref   = [pi; pi; 0; 0];

% ---- tracking error ----
    e    = x - xref;
    e(1) = mod(e(1)+pi,2*pi)-pi;
    e(2) = mod(e(2)+pi,2*pi)-pi;

% ---- integral sliding variable ----                                  NEW
    s    = c.'*e + ki*xi;          % was: s = c.'*e

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

% ---- SMC law (equivalent + switching) ----
    cg  = (pu/detM)*( c(3) - c(4)*M12 );          % c'*g(x)
    ueq = -(c.'*f0 + ki*e(1))/cg;                 % +ki*e1 from d/dt of xi    NEW
    usw = -(k/cg)*max(min(s/phi,1),-1);
    u   = max(min(ueq + usw, 1), -1);             % scalar, saturated

% ---- conditional integration ----
    if abs(ueq + usw) < 1                          % integrate only if not saturating
        xi = xi + h*e(1);                          % forward-Euler integral of e1
        xi = max(min(xi, xi_max), -xi_max);        % clamp
    end
end