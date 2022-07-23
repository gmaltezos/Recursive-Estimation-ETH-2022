function [postParticles] = Estimator(prevPostParticles, sens, act, estConst, km)
% The estimator function. The function will be called in two different
% modes: If km==0, the estimator is initialized. If km > 0, the
% estimator does an iteration for a single sample time interval using the 
% previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurement and control inputs.
%
% Inputs:
%   prevPostParticles   previous posterior particles at time step k-1
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .rho: wall p1 to p2 offset [m]
%                       .kappa: wall p8 to p9 offset [m]
%                           
%   sens                Sensor measurement z(k), scalar
%
%   act                 Control inputs u(k-1), [1x2]-vector
%                       act(1): u_f, forward control input
%                       act(2): u_phi, angular control input
%
%   estConst            estimator constants (as in EstimatorConst.m)
%
%   km                  time index k, scalar
%                       corresponds to continous time t = k*Ts
%                       If km==0 initialization, otherwise estimator
%                       iteration step.
%
% Outputs:
%   postParticles       Posterior particles at time step k
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .rho: wall p1 to p2 offset [m]
%                       .kappa: wall p8 to p9 offset [m]
%
%
% Class:
% Recursive Estimation
% Spring 2022
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Enrico Mion, Bhavya Sukhija, Jin Cheng
% enmion@ethz.ch
% sukhijab@ethz.ch
% jicheng@student.ethz.ch

% Set number of particles:
N_particles = 3500; % obviously, you will need more particles than 10.

%% Mode 1: Initialization
if (km == 0)
    % Do the initialization of your estimator here!
    % Distribution Calculations
    phiA = 2*pi*rand(1,N_particles/2);
    phiB = 2*pi*rand(1,N_particles/2);
    rA = estConst.d*sqrt(rand(1,N_particles/2));
    rB = estConst.d*sqrt(rand(1,N_particles/2));
    postParticles.x_r = [estConst.pA(1) + rA.*cos(phiA), estConst.pB(1) + rB.*cos(phiB)]; % 1xN_particles matrix
    postParticles.y_r = [estConst.pA(2) + rA.*sin(phiA), estConst.pB(2) + rB.*sin(phiB)]; % 1xN_particles matrix

    postParticles.phi = -estConst.phi_0 + 2*estConst.phi_0*rand(1,N_particles); % 1xN_particles matrix uniform distribution with mean 0
    postParticles.rho = -estConst.m + 2*estConst.m*rand(1,N_particles); % 1xN_particles matrix uniform distribution with mean 0
    postParticles.kappa = -estConst.l + 2*estConst.l*rand(1,N_particles); % 1xN_particles matrix uniform distribution with mean 0
    
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.

% Implement your estimator here!
% Prior Update:
% Inputs and constants
uf = act(1) * ones(1,N_particles);
uphi = act(2) * ones(1,N_particles);
sigmaf = estConst.sigma_f/2;
sigmaphi = estConst.sigma_phi/2;

% Noise Samples
vf = -sigmaf + 2*sigmaf*rand(1,N_particles);
vphi = -sigmaphi + 2*sigmaphi*rand(1,N_particles);

% Apply the process Equation to the particles
x_pp = prevPostParticles.x_r + (uf + vf).*cos(prevPostParticles.phi);
y_pp = prevPostParticles.y_r + (uf + vf).*sin(prevPostParticles.phi);
phi_pp = prevPostParticles.phi + uphi + vphi;
rho_pp = prevPostParticles.rho;
k_pp = prevPostParticles.kappa;

% Posterior Update:
% Constants
epsilon = estConst.epsilon;
xCont = estConst.contour(:,1);
yCont = estConst.contour(:,2);

w = zeros(1,N_particles);

% Contour Lines (By using consecutive points its time)
xContall = repmat(xCont,1,N_particles);
xContall(8,:) =  k_pp;
xContall(9,:) = k_pp;
yContall = repmat(yCont,1,N_particles);
yContall(1,:) =  rho_pp;
yContall(2,:) = rho_pp;
xContall2 = [xContall(2:10,:); xContall(1,:)];
yContall2 = [yContall(2:10,:); yContall(1,:)];

% Robot Line Direction
aRob = tan(phi_pp);
bRob = y_pp - aRob.*x_pp;

% Checking the interesction of the line of the robot with the line
% segments of the contour
t = (bRob - yContall + aRob.*xContall)./(yContall2 - yContall + aRob.*(xContall - xContall2));
for i=1:N_particles
    xinter = xContall(t(:,i)>=0 & t(:,i)<=1,i) + (xContall2(t(:,i)>=0 & t(:,i)<=1,i)- xContall(t(:,i)>=0 & t(:,i)<=1,i)).*t(t(:,i)>=0 & t(:,i)<=1,i);
    yinter = yContall(t(:,i)>=0 & t(:,i)<=1,i) +(yContall2(t(:,i)>=0 & t(:,i)<=1,i) - yContall(t(:,i)>=0 & t(:,i)<=1,i)).*t(t(:,i)>=0 & t(:,i)<=1,i);

    % Checking if the points lie in front of the robot
    if 0 <= phi_pp(i) && phi_pp(i) < pi/2
        xinter_kept = xinter(xinter > x_pp(i) & yinter > y_pp(i));
        yinter_kept = yinter(xinter > x_pp(i) & yinter > y_pp(i));
    elseif pi/2 <= phi_pp(i) && phi_pp(i) < pi
        xinter_kept = xinter(xinter < x_pp(i) & yinter > y_pp(i));
        yinter_kept = yinter(xinter < x_pp(i) & yinter > y_pp(i));
    elseif pi <= phi_pp(i) && phi_pp(i) < 3*pi/2
        xinter_kept = xinter(xinter < x_pp(i) & yinter < y_pp(i));
        yinter_kept = yinter(xinter < x_pp(i) & yinter < y_pp(i));  
    else
        xinter_kept = xinter(xinter > x_pp(i) & yinter < y_pp(i));
        yinter_kept = yinter(xinter > x_pp(i) & yinter < y_pp(i));
    end

    % Finding the final closest point
    dist = sqrt((yinter_kept - y_pp(i)).^2 + (xinter_kept - x_pp(i)).^2);
    [mindist,idx] = min(dist); 
%     xinter_final = xinter_kept(idx);
%     yinter_final = yinter_kept(idx);
    
    if isempty(idx)
        w(i) = 100; % Something big, so this particle to not be considered
    else
        % Draw samples from measurement noise with the measurement equation
        w(i) = abs(sens - mindist); % Symmetric pdf
    end
end

% Calculating Measurement likelihood
p_w = zeros(size(w));
if epsilon == 0
    p_w(w == 0) = 1000; % High probability of w being zero
else
    p_w(w >= 0 & w < 2*epsilon) = -w(w >= 0 & w < 2*epsilon)/(5*(epsilon^2)) + 2/(5*epsilon);       
    p_w(w >= 2*epsilon & w < 2.5*epsilon) = w(w >= 2*epsilon & w < 2.5*epsilon)/(2.5*(epsilon^2)) - 2/(2.5*epsilon);
    p_w(w >= 2.5*epsilon & w < 3*epsilon) = -w(w >= 2.5*epsilon & w < 3*epsilon)/(2.5*(epsilon^2)) + 3/(2.5*epsilon);
end

if(sum(p_w) == 0)
%     warning('Particle weights were all zero')

    postParticles.x_r = x_pp;
    postParticles.y_r = y_pp;
    postParticles.phi = phi_pp;
    postParticles.rho = rho_pp;
    postParticles.kappa = k_pp;
else
    betas = p_w/sum(p_w);

    % Resampling
    betaCumSum = cumsum(betas);
    randNumber = rand(1,N_particles);
    for i = 1:N_particles
        ind = find(betaCumSum >= randNumber(i),1,'first');
        postParticles.x_r(i) = x_pp(ind);
        postParticles.y_r(i) = y_pp(ind);
        postParticles.phi(i) = phi_pp(ind);
        postParticles.rho(i) = rho_pp(ind);
        postParticles.kappa(i) = k_pp(ind);
    end
end

% Roughening
K = 0.02; % Roughening parameter, tuned for the specific number of particles used
D = 5;    % Dimension of the state space
E = [max(postParticles.x_r) - min(postParticles.x_r);
     max(postParticles.y_r) - min(postParticles.y_r);
     max(postParticles.phi) - min(postParticles.phi);
     max(postParticles.rho) - min(postParticles.rho);
     max(postParticles.kappa) - min(postParticles.kappa)];
 StdRough = K.*diag(E)*N_particles^(-1/D);
 deltaX = StdRough*randn(D,N_particles);

postParticles.x_r = postParticles.x_r + deltaX(1,:);
postParticles.y_r = postParticles.y_r + deltaX(2,:);
postParticles.phi = postParticles.phi + deltaX(3,:);
postParticles.rho = postParticles.rho + deltaX(4,:);
postParticles.kappa = postParticles.kappa + deltaX(5,:);

end % end estimator