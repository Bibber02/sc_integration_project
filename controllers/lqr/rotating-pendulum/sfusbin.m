function [sys,x0,str,ts] = sfusbin(t,x,u,flag,Ts,h,channels)
ensureRotpendHardwarePath(mfilename('fullpath'));
[sys,x0,str,ts] = rotpend_sfusbin(t,x,u,flag,Ts,h,channels);
end
