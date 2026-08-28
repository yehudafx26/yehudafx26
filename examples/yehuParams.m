function [fStringParams, cStringParams, bodyParams, fingerParams, bowParams, frictionParams] = yehuParams()
% YEHUPARAMS.M - Provides physical parameters for yehu simulation 
%
% STRING PARAMETERS
% [f0, L, r, rho, E, sig0, sig1, xcoup]
% fundamental frequency [Hz], length [m], density [kg/m^3], 
% Young's modulus [Pa], freq. independant damping [1/s], 
% freq. dependant damping [m^2/s], normalized bridge position [-]
%
% FINGER PARAMETERS
% [mf, kf, rf, alf, uf0]
% mass [kg], stiffness [N/m^alf], damping [s/m], initial position [m]
%
% BOW PARAMETERS
% [mh, kh, rh, xb]
% bow hair mass [kg], stiffness [N/m], damping [kg/s], 
% normalized bow position [-]
%
% FRICTION PARAMETERS
% [s0, s1_bar, vS, muC, muS, Sexp]
% bristle stiffness [N/m], bristle damping [kg/s], Stribeck velocity,
% Columb friction coefficient, Stribeck friction ceofficient, Stribeck
% exponent
%
% BODY PARAMETERS
% [Mq, Kq, Rq]
% Modal mass [kg], modal stiffness [N/m], modal dampings [1/s]

%% STRING PARAMETERS
% [f0, L, r, rho, E, sig0, sig1]
L = 0.38;
xcoup = (L - 0.08) / L; % bridge coupling position 
fStringParams = [ 344, L, 5.5e-4, 1300, 10e9, 0.9754, 0.0026, xcoup];
cStringParams = [ 513, L, 4e-4, 1300, 10e9, 1.5884, 0.0032, xcoup ];


%% FINGER PARAMETERS
% [mf, kf, rf, alf, uf0]
uf0 = 5e-4; % initial hammer position 
fingerParams = [0.02, 5e5, 20, 2.3, uf0];

%% BOW PARAMETERS
% [mh, kh, rh, xb]
xb = 0.255 / L ; % bowing position 
bowParams = [4.5e-3, 4.8297e4, 5.7674, xb];

%% FRICTION PARAMETERS
% [s0, s1_bar, vS, muC, muS, Sexp]
frictionParams = [3.1860e4, 0.0027, 0.2280, 0.5071, 1.0207, 2];

%% BODY PARAMETERS
% [Mq, Kq, Rq]
load("yehuBridgeAdmittance.mat");
bodyParams = [M, K, R];
end