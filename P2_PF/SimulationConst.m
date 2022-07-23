function const = SimulationConst()
% 
% Define the constants used in the simulation.  These constants are not 
% accessible to the estimator.
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

%% The wall contour - (x,y) coordinates of corners as in Table 1
const.contour = [0.50, 0.00;
                 2.50, 0.00;
                 2.50, 1.50;
                 3.00, 2.00;
                 2.00, 3.00;
                 1.25, 2.25;
                 1.00, 2.50;
                 0.00, 2.50;
                 0.00, 0.50;
                 0.50, 0.50];
             
%% Initialization
const.pA = [1.4,0.8]; % Center point pA of the initial position distribution
const.pB = [2.0,1.8]; % Center point pB of the initial position distribution
const.d = 0.1;  % Radius of shaded regions for initialization

const.phi_0 = pi/3; % Initial heading is uniformly distr. in [-phi_0,phi_0]

const.m = 0.1;  % Uniform distribution in [0.1,0.2] parameter of points p1 and p2
const.l = 0.2;  % Uniform distribution in [0.1,0.2] parameter of points p8 and p9

%% Noise properties
% process noise
const.sigma_phi = 0.05; % Parameter for process noise v_phi
const.sigma_f = 0.01; % Parameter for process noise v_f

% measurement noise
const.epsilon = 0.01; % Parameter for measurement noise w

%% Times
% Number of samples (discrete time steps) of the simulation.
const.N = 500;
