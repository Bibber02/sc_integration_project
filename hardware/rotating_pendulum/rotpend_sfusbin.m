function [sys,x0,str,ts] = rotpend_sfusbin(t,x,u,flag,Ts,h,channels)
%ROTPEND_SFUSBIN USB input from a DCSC setup into Simulink.

switch flag
    case 0
        [sys,x0,str,ts] = mdlInitializeSizes(Ts,h);
    case 3
        sys = mdlOutputs(u,x,h,channels);
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
sizes.NumOutputs     = 7;
sizes.NumInputs      = 0;
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

function sys = mdlOutputs(u,x,h,channels) %#ok<INUSD>
d = fugiboard('Read', h);
s = rem(floor(d(8) * [0.5 0.25 0.125 0.0625]), 2)';
sys = [d(channels); s; d(4)];

function sys = mdlTerminate(h)
fugiboard('Write', h, 0, 1, 0, 0);
fugiboard('Close', h);
sys = [];
