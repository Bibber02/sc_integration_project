function [sys,x0,str,ts] = rotpend_sfusbout(t,x,u,flag,Ts,h)
%ROTPEND_SFUSBOUT USB output from Simulink to a DCSC setup.

switch flag
    case 0
        [sys,x0,str,ts] = mdlInitializeSizes(Ts,h);
    case 3
        sys = mdlOutputs(u,x,Ts,h);
        x0 = [];
        str = [];
        ts = [];
    case 9
        sys = mdlTerminate(h);
        x0 = [];
        str = [];
        ts = [];
    otherwise
        sys = [];
        x0 = [];
        str = [];
        ts = [];
end

function [sys,x0,str,ts] = mdlInitializeSizes(Ts,h) %#ok<INUSD>
sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 0;
sizes.NumInputs      = 1;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0 = [];
str = [];
ts = [-1 0];
if length(Ts) > 0
    ts(1) = Ts(1);
end
if length(Ts) > 1
    ts(2) = Ts(2);
end
tic;

function sys = mdlOutputs(u,x,Ts,h) %#ok<INUSD>
fugiboard('Write', h, 0, 1, u, 0);
while toc < Ts
end
tic;
sys = [];

function sys = mdlTerminate(h)
fugiboard('Write', h, 0, 1, 0, 0);
fugiboard('Close', h);
sys = [];
