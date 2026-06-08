function fugihandle = rotpend_calib()
%ROTPEND_CALIB Initialize the FPGA interface board.

fugihandle = fugiboard('Open', 'Pendulum1');
fugihandle.WatchdogTimeout = 0.5;
fugiboard('SetParams', fugihandle);
fugiboard('Write', fugihandle, 0, 0, 0, 0);
fugiboard('Write', fugihandle, 5, 1, 0, 0);
data = fugiboard('Read', fugihandle);
model = bitshift(data(1), -5);
version = bitand(data(1), 31);
disp(sprintf('FPGA setup %d, version %d', model, version));
fugiboard('Write', fugihandle, 0, 1, 0, 0);
pause(0.1);
end
