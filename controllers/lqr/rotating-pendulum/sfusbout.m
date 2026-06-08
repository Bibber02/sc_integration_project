function [sys,x0,str,ts] = sfusbout(t,x,u,flag,Ts,h)
ensureRotpendHardwarePath(mfilename('fullpath'));
[sys,x0,str,ts] = rotpend_sfusbout(t,x,u,flag,Ts,h);
end
