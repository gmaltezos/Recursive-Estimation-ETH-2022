function [posEst,linVelEst,oriEst,windEst,driftEst,...
    posVar,linVelVar,oriVar,windVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,windEst,driftEst,...
%    posVar,linVelVar,oriVar,windVar,driftVar,estState] =
% Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time t_k, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConst.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   windEst         wind direction estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   windVar         variance of wind direction estimate(time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2022
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Enrico Mion, Bhavya Sukhija, Jin Cheng
% enmion@ethz.ch
% sukhijab@ethz.ch
% jicheng@student.ethz.ch

%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!

    % initial state mean
    posEst = [0, 0]; % 1x2 matrix
    linVelEst = [0, 0]; % 1x2 matrix
    oriEst = 0; % 1x1 matrix
    windEst = 0; % 1x1 matrix
    driftEst = 0; % 1x1 matrix

    % initial state variance
    posVar = [estConst.StartRadiusBound^2/4, estConst.StartRadiusBound^2/4]; % 1x2 matrix
    linVelVar = [0, 0]; % 1x2 matrix known initial velocity
    oriVar = (estConst.RotationStartBound+estConst.RotationStartBound)^2/12; % 1x1 matrix from uniform distribution
    windVar = (estConst.WindAngleStartBound+estConst.WindAngleStartBound)^2/12; % 1x1 matrix from uniform distribution
    driftVar = 0;

    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar, linVelVar, oriVar, windVar, driftVar]);
    % estimator state
    estState.xm = [posEst, linVelEst, oriEst, windEst, driftEst]';
    % time of last update
    estState.tm = tm;
    return;
end
% Estimator iteration.
% get time since last estimator update
dt = tm - estState.tm;
estState.tm = tm; % update measurement update time

% prior update
% Constants
Cr = estConst.rudderCoefficient;
ut = actuate(1);
ur = actuate(2);
Cda = estConst.dragCoefficientAir;
Cdh = estConst.dragCoefficientHydr;
Cw = estConst.windVel;
Q = diag([estConst.DragNoise, estConst.RudderNoise, estConst.WindAngleNoise, estConst.GyroDriftNoise]);

%%% Symbolic calculations for A and L matrices %%%
%     syms px1 py1 sx1 sy1 fi1 ro1 bi1 Cr ut ur Cda Cdh Cw vd vr vw vb
%     f1 = sx1;
%     f2 = sy1;
%     f3 = cos(fi1)*(tanh(ut)-Cdh*(sx1^2+sy1^2)*(1+vd)) - Cda*(sx1-Cw*cos(ro1))*sqrt((sx1 - Cw*cos(ro1))^2+(sy1-Cw*sin(ro1))^2);
%     f4= sin(fi1)*(tanh(ut)-Cdh*(sx1^2+sy1^2)*(1+vd)) - Cda*(sy1-Cw*sin(ro1))*sqrt((sx1 - Cw*cos(ro1))^2+(sy1-Cw*sin(ro1))^2);
%     f5=Cr*ur*(1+vr);
%     f6=vw;
%     f7=vb;
%     jacobian([f1,f2,f3,f4,f5,f6,f7],[px1,py1,sx1,sy1,fi1,ro1,bi1])
%     jacobian([f1,f2,f3,f4,f5,f6,f7],[vd, vr, vw, vb])

% Solve the dynamics
q = @my_ode_dynamics;
[~,xp] = ode45(q, [tm-dt tm], estState.xm);
xp = xp(end,:)';

% Calculating the linearized dynamics
A = zeros(7,7);
A(1,3) = 1;
A(2,4) = 1;
A(3,3) = - Cda*((xp(3) - Cw*cos(xp(6)))^2 + (xp(4) - Cw*sin(xp(6)))^2)^(1/2) - 2*Cdh*xp(3)*cos(xp(5)) - (Cda*(xp(3) - Cw*cos(xp(6)))*(2*xp(3) - 2*Cw*cos(xp(6))))/(2*((xp(3) - Cw*cos(xp(6)))^2 + (xp(4) - Cw*sin(xp(6)))^2)^(1/2));
A(4,3) = - 2*Cdh*xp(3)*sin(xp(5)) - (Cda*(xp(4) - Cw*sin(xp(6)))*(2*xp(3) - 2*Cw*cos(xp(6))))/(2*((xp(3) - Cw*cos(xp(6)))^2 + (xp(4) - Cw*sin(xp(6)))^2)^(1/2));
A(3,4) = - 2*Cdh*xp(4)*cos(xp(5)) - (Cda*(xp(3) - Cw*cos(xp(6)))*(2*xp(4) - 2*Cw*sin(xp(6))))/(2*((xp(3) - Cw*cos(xp(6)))^2 + (xp(4) - Cw*sin(xp(6)))^2)^(1/2));
A(4,4) = - Cda*((xp(3) - Cw*cos(xp(6)))^2 + (xp(4) - Cw*sin(xp(6)))^2)^(1/2) - 2*Cdh*xp(4)*sin(xp(5)) - (Cda*(xp(4) - Cw*sin(xp(6)))*(2*xp(4) - 2*Cw*sin(xp(6))))/(2*((xp(3) - Cw*cos(xp(6)))^2 + (xp(4) - Cw*sin(xp(6)))^2)^(1/2));
A(3,5) = -sin(xp(5))*(tanh(ut) - Cdh*(xp(3)^2 + xp(4)^2));
A(4,5) = cos(xp(5))*(tanh(ut) - Cdh*(xp(3)^2 + xp(4)^2));
A(3,6) = - (Cda*(2*Cw*sin(xp(6))*(xp(3) - Cw*cos(xp(6))) - 2*Cw*cos(xp(6))*(xp(4) - Cw*sin(xp(6))))*(xp(3) - Cw*cos(xp(6))))/(2*((xp(3) - Cw*cos(xp(6)))^2 + (xp(4) - Cw*sin(xp(6)))^2)^(1/2)) - Cda*Cw*sin(xp(6))*((xp(3) - Cw*cos(xp(6)))^2 + (xp(4) - Cw*sin(xp(6)))^2)^(1/2);
A(4,6) = Cda*Cw*cos(xp(6))*((xp(3) - Cw*cos(xp(6)))^2 + (xp(4) - Cw*sin(xp(6)))^2)^(1/2) - (Cda*(2*Cw*sin(xp(6))*(xp(3) - Cw*cos(xp(6))) - 2*Cw*cos(xp(6))*(xp(4) - Cw*sin(xp(6))))*(xp(4) - Cw*sin(xp(6))))/(2*((xp(3) - Cw*cos(xp(6)))^2 + (xp(4) - Cw*sin(xp(6)))^2)^(1/2));

L = zeros(7,4);
L(3,1) = -Cdh*cos(xp(5))*(xp(3)^2 + xp(4)^2);
L(4,1) = -Cdh*sin(xp(5))*(xp(3)^2 + xp(4)^2);
L(5,2) = Cr*ur;
L(6,3) = 1;
L(7,4) = 1;

% Solving for the variance
Ppeq = @(t,P) my_ode_var(t,P,A,L,Q);
[~,Pp] = ode45(Ppeq,[tm-dt tm],estState.Pm(:));
Pp = reshape(Pp(end,:),7,7,1);

%%% Measurement Update
% Constants
xa = estConst.pos_radioA(1);
ya = estConst.pos_radioA(2);
xb = estConst.pos_radioB(1);
yb = estConst.pos_radioB(2);
xc = estConst.pos_radioC(1);
yc = estConst.pos_radioC(2);
R = diag([estConst.DistNoiseA, estConst.DistNoiseB, estConst.DistNoiseC, estConst.GyroNoise, estConst.CompassNoise]);

% Measurement Dynamics
%%%% Symbolic calculations for the measurement update
%     syms px1 py1 sx1 sy1 fi1 ro1 bi1 Cr ut ur Cda Cdh Cw vd vr vw vb xa xb xc ya yb yc
%     f1 = ((px1 - xa)^2 + (py1 - ya)^2)^(1/2);
%     f2 = ((px1 - xb)^2 + (py1 - yb)^2)^(1/2);
%     f3 = ((px1 - xc)^2 + (py1 - yc)^2)^(1/2);
%     f4 = fi1 + bi1;
%     f5 = fi1;
%     jacobian([f1,f2,f3,f4,f5],[px1,py1,sx1,sy1,fi1,ro1,bi1])
%     jacobian([f1,f2,f3,f4,f5],[vd, vr, vw, vb])

% Measurement Equations
za = ((xp(1) - xa)^2 + (xp(2) - ya)^2)^(1/2);
zb = ((xp(1) - xb)^2 + (xp(2) - yb)^2)^(1/2);
zc = ((xp(1) - xc)^2 + (xp(2) - yc)^2)^(1/2);
zg = xp(5) + xp(7);
zn = xp(5);
h = [za, zb, zc, zg, zn]';

% Linearized Measurement Dynamics
H = zeros(5,7);
H(1,1) = (2*xp(1) - 2*xa)/(2*((xp(1) - xa)^2 + (xp(2) - ya)^2)^(1/2));
H(1,2) = (2*xp(2) - 2*ya)/(2*((xp(1) - xa)^2 + (xp(2) - ya)^2)^(1/2));
H(2,1) = (2*xp(1) - 2*xb)/(2*((xp(1) - xb)^2 + (xp(2) - yb)^2)^(1/2));
H(2,2) = (2*xp(2) - 2*yb)/(2*((xp(1) - xb)^2 + (xp(2) - yb)^2)^(1/2));
H(3,1) = (2*xp(1) - 2*xc)/(2*((xp(1) - xc)^2 + (xp(2) - yc)^2)^(1/2));
H(3,2) = (2*xp(2) - 2*yc)/(2*((xp(1) - xc)^2 + (xp(2) - yc)^2)^(1/2));
H(4,5) = 1;
H(4,7) = 1;
H(5,5) = 1;

M = eye(5,5);

% ifclauses to check if you have measurements in this particular step
if isinf(sense(1))
    sense(1) = [];
    R(1,:) = [];
    R(:,1) = [];
    h(1) = [];
    H(1,:) = [];
    M(1,:) = [];
    M(:,1) = [];
end

if length(sense) < 5 && isinf(sense(1))
    sense(1) = [];
    R(1,:) = [];
    R(:,1) = [];
    h(1) = [];
    H(1,:) = [];
    M(1,:) = [];
    M(:,1) = [];
elseif length(sense) >= 5 && isinf(sense(2))
    sense(2) = [];
    R(2,:) = [];
    R(:,2) = [];
    h(2) = [];
    H(2,:) = [];
    M(2,:) = [];
    M(:,2) = [];
end

if length(sense) < 4 && isinf(sense(1))
    sense(1) = [];
    R(1,:) = [];
    R(:,1) = [];
    h(1) = [];
    H(1,:) = [];
    M(1,:) = [];
    M(:,1) = [];
elseif length(sense) < 5 && isinf(sense(2))
    sense(2) = [];
    R(2,:) = [];
    R(:,2) = [];
    h(2) = [];
    H(2,:) = [];
    M(2,:) = [];
    M(:,2) = [];
elseif length(sense) >= 5 && isinf(sense(3))
    sense(3) = [];
    R(3,:) = [];
    R(:,3) = [];
    h(3) = [];
    H(3,:) = [];
    M(3,:) = [];
    M(:,3) = [];
end

% Kalman Filter Update
K = Pp*H'/(H*Pp*H'+M*R*M');
estState.xm = xp + K*(sense' - h);
estState.Pm = (eye(size(K,1)) - K*H)*Pp;

% Get resulting estimates and variances
% Output quantities
posEst = estState.xm(1:2);
linVelEst = estState.xm(3:4);
oriEst = estState.xm(5);
windEst = estState.xm(6);
driftEst = estState.xm(7);

posVar = [estState.Pm(1,1), estState.Pm(2,2)];
linVelVar = [estState.Pm(3,3), estState.Pm(4,4)];
oriVar = estState.Pm(5,5);
windVar = estState.Pm(6,6);
driftVar = estState.Pm(7,7);

% end
% Definition of functions for solving with ODE
    function q = my_ode_dynamics(t,x)
        pxdot = x(3);
        pydot = x(4);
        sxdot = cos(x(5))*(tanh(ut)-Cdh*(x(3)^2+x(4)^2)) - Cda*(x(3)-Cw*cos(x(6)))*sqrt((x(3) - Cw*cos(x(6)))^2+(x(4)-Cw*sin(x(6)))^2);
        sydot = sin(x(5))*(tanh(ut)-Cdh*(x(3)^2+x(4)^2)) - Cda*(x(4)-Cw*sin(x(6)))*sqrt((x(3) - Cw*cos(x(6)))^2+(x(4)-Cw*sin(x(6)))^2);
        q = [pxdot; pydot; sxdot; sydot; Cr*ur; 0; 0];
    end

    function Ppeq = my_ode_var(t,P,A,L,Q)
        P = reshape(P,size(A));
        Ppeq = A*P + P*A' + L*Q*L';
        Ppeq = Ppeq(:);
    end
end
