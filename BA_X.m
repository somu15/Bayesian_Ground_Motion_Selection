%% Sets up the predictor variable matrix for computing Bayesian CS
%  (GMPM is Boore-Atkinson 2008; small and large magnitude earthquakes)

% "A Bayesian Treatment of the Conditional Spectrum Approach for Ground
% Motion Selection". Report by Somayajulu Dhulipala, Jack Baker and
% Madeleine Flint.
%
% Created by Somayajulu Dhulipala, 20 May 2017
% Modified by Somayajulu Dhulipala, 23 December 2017

%% INPUTS

% magnitude      =      the earthquake magnitude vector (N X 1)
% Rjb            =      the Joyner-Boore distance (Km) vector (N X 1)
% Vs             =      the shear-wave velocity vector (m/s) (N X 1)
% SS             =      indicator variable for Strike-Slip fault (N X 1)
% N              =      indicator variable for Normal fault (N X 1)
% R              =      indicator variable for Reverse fault (N X 1)
% U              =      indicator variable for Unspecified fault (N X 1)
% M_ref          =      reference magnitude (1 X 1)
% R_ref          =      reference distance (1 X 1)
% M_h            =      hinge magnitude (1 X 1)
% Vs_ref         =      reference shear-wave velocity (1 X 1)
% h              =      Boore-Atkinson 2008 GMPM variable (1 X 1)

%% OUTPUTS

% X              =      predictor variable matrix for Boore-Atkinson 2008
% (N X 11)

%% EXPLICITLY DEFINED FUNCTIONS/ROUTINES USED

% NONE

%% FUNCTION BEGINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X] = BA_X(magnitude, Rjb, Vs, SS, N, R, U, M_ref, R_ref, M_h, Vs_ref, h)


index = find(magnitude <= M_h);
e_5 = zeros(length(magnitude),1);
e_5(index) = 1;
e_5 = e_5.*(magnitude-M_h);
e_6 = zeros(length(magnitude),1);
e_6(index) = 1;
e_6 = e_6.*(magnitude-M_h).^2;


index = find(magnitude > M_h);
e_7 = zeros(length(magnitude),1);
e_7(index) = 1;
e_7 = e_7.*(magnitude-M_h);

R_R = sqrt(Rjb.^2+h^2);

X = [U SS N R e_5 e_6 e_7 log(R_R/R_ref)...
    (magnitude-M_ref).*log(R_R/R_ref) R_R-R_ref log(Vs/Vs_ref)];


end