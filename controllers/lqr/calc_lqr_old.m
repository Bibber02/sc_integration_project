function [K_lqr, closedLoopPoles] = calc_lqr(sys_disc, Q_lqr, R_lqr)
    % Extract discrete matrices from the LTI system
    Ad = sys_disc.A;
    Bd = sys_disc.B;

    % Discrete-time LQR gain for model coordinates:
    %   u_model_dev = -K_lqr * x_dev.
    % The real hardware command may need an extra sign conversion; setup_lqr.m
    % computes K_lqr_command for the Simulink Gain block.
    [K_lqr, ~, closedLoopPoles] = dlqr(Ad, Bd, Q_lqr, R_lqr);
    disp('Closed-loop poles:');
    disp(closedLoopPoles);

    disp('K_lqr = ');
    disp(K_lqr);

end
