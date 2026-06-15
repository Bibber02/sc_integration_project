function [K_lqr, closedLoopPoles, lqrInfo] = calc_lqr(sys_disc, Q_lqr, R_lqr, Ts)
    %CALC_LQR Discrete-time LQR design with useful diagnostics.
    %
    % Model convention:
    %   x_dev(k+1)  = Ad*x_dev(k) + Bd*u_model_dev(k)
    %   u_model_dev = -K_lqr*x_dev
    %
    % If the real command input has the opposite sign, do not change the
    % LQR design here. Instead convert the sign in setup_lqr.m.

    if nargin < 4 || isempty(Ts)
        Ts = sys_disc.Ts;
    end

    % Extract discrete matrices from the LTI system.
    Ad = sys_disc.A;
    Bd = sys_disc.B;

    n = size(Ad,1);

    % Discrete-time LQR gain.
    [K_lqr, S_lqr, closedLoopPoles] = dlqr(Ad, Bd, Q_lqr, R_lqr);

    Acl = Ad - Bd*K_lqr;
    openLoopPoles = eig(Ad);
    closedLoopPoles = eig(Acl);
    maxClosedLoopPoleMagnitude = max(abs(closedLoopPoles));
    rankCtrb = rank(ctrb(Ad,Bd));

    if ~isempty(Ts) && Ts > 0
        s_closedLoopPoles = log(closedLoopPoles)/Ts;
    else
        s_closedLoopPoles = NaN(size(closedLoopPoles));
    end

    lqrInfo = struct();
    lqrInfo.Ad = Ad;
    lqrInfo.Bd = Bd;
    lqrInfo.Acl = Acl;
    lqrInfo.S_lqr = S_lqr;
    lqrInfo.openLoopPoles = openLoopPoles;
    lqrInfo.closedLoopPoles = closedLoopPoles;
    lqrInfo.s_closedLoopPoles = s_closedLoopPoles;
    lqrInfo.maxClosedLoopPoleMagnitude = maxClosedLoopPoleMagnitude;
    lqrInfo.rankCtrb = rankCtrb;
    lqrInfo.nStates = n;

    fprintf('\ncalc_lqr diagnostics:\n');
    fprintf('  Controllability rank:        %d / %d\n', rankCtrb, n);
    fprintf('  Max closed-loop pole |z|:    %.6f\n', maxClosedLoopPoleMagnitude);

    disp('  Open-loop poles:');
    disp(openLoopPoles);

    disp('  Closed-loop poles:');
    disp(closedLoopPoles);

    disp('  K_lqr = ');
    disp(K_lqr);
end
