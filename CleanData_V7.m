%% Cleans the ground motion database

% "A Bayesian Treatment of the Conditional Spectrum Approach for Ground
% Motion Selection". Report by Somayajulu Dhulipala, Jack Baker and
% Madeleine Flint.
%
% Created by Somayajulu Dhulipala, 20 May 2017
% Modified by Somayajulu Dhulipala, 23 December 2017

%% INPUTS

% IM_matrix      =      the intensity measure matrix (N X (Nt+Nn))
% N is the number of observations, Nt is the number of time periods
% Nn is the number of non-spectral IMs
% magnitude      =      the earthquake magnitude vector (N X 1)
% Rjb            =      the Joyner-Boore distance (Km) vector (N X 1)
% Vs             =      shear-wave velocity (m/s) vector (N X 1)
% mechanism      =      fault-mechanism vector (N X 1)
% req_per        =      the required time periods vector (Nt X 1)
% restr_vec      =      non-spectral IM level restriction matrix (Nn X 2)
% for example, if there are two non-spectral IMs PGA and PGV and their
% respective minimum and maximum values are (0.01 2g) and (1 200cm/s),
% then restr_vec is [0.01 2;1 200]

% PEER fault mechanism description
% 0 - Strike slip
% 1 - Normal slip
% 2 - Reverse slip
% 3 - Reverse oblique
% 4 - Normal oblique

%% OUTPUTS

% IM_matrix      =      the cleaned intensity measure matrix (N X (Nt+Nn))
% magnitude      =      the cleaned earthquake magnitude vector (N X 1)
% Rjb            =      cleaned Joyner-Boore distance (Km) vector (N X 1)
% Vs             =      cleaned shear-wave velocity (m/s) vector (N X 1)
% mechanism      =      cleaned fault-mechanism vector (N X 1)
% SS             =      Strike-Slip fault indicator vector (N X 1)
% N              =      Normal fault indicator vector (N X 1)
% R              =      Reverse fault indicator vector (N X 1)
% U              =      Unspecified fault indicator vector (N X 1)
% mechanism_BA   =      Boore-Atkinson type mechanism indicator (NA)

%% EXPLICITLY DEFINED FUNCTIONS/ROUTINES USED

% NONE

%% FUNCTION BEGINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [IM_matrix, magnitude, Rjb, Vs, mechanism, SS, N, R, U, mechanism_BA] = CleanData_V7(IM_matrix, magnitude, Rjb, Vs, mechanism, req_per, restr_vec)

siz = size(IM_matrix);

% IM restrictions

index = find(IM_matrix(:,1) > 2);
IM_matrix(index,:) = [];
magnitude(index) = [];
Rjb(index,:) = [];
Vs(index,:) = [];
mechanism(index,:) = [];

index = find(IM_matrix(:,17) < 0.01);
IM_matrix(index,:) = [];
magnitude(index) = [];
Rjb(index,:) = [];
Vs(index,:) = [];
mechanism(index,:) = [];

if siz(2) > length(req_per)
    for ii = (length(req_per)+1):siz(2)
        
        temp_ind = length(req_per)+(ii-length(req_per));
index = find(IM_matrix(:, temp_ind) > restr_vec((ii-length(req_per)),2));
IM_matrix(index,:) = [];
magnitude(index) = [];
Rjb(index,:) = [];
Vs(index,:) = [];
mechanism(index,:) = [];

index = find(IM_matrix(:, temp_ind) < restr_vec((ii-length(req_per)),1));
IM_matrix(index,:) = [];
magnitude(index) = [];
Rjb(index,:) = [];
Vs(index,:) = [];
mechanism(index,:) = [];

    end
end

% Magnitude restrictions

index = find(magnitude < 4);
IM_matrix(index,:) = [];
magnitude(index) = [];
Rjb(index,:) = [];
Vs(index,:) = [];
mechanism(index,:) = [];

% Distance restrictions

index = find(Rjb > 200);
IM_matrix(index,:) = [];
magnitude(index) = [];
Rjb(index,:) = [];
Vs(index,:) = [];
mechanism(index,:) = [];

index = find(Rjb < 1);
IM_matrix(index,:) = [];
magnitude(index) = [];
Rjb(index,:) = [];
Vs(index,:) = [];
mechanism(index,:) = [];

% Vs30 restrictions

index = find(Vs < 0);
IM_matrix(index,:) = [];
magnitude(index) = [];
Rjb(index,:) = [];
Vs(index,:) = [];
mechanism(index,:) = [];

% Mechanism restrictions

index = find(mechanism < 0);
IM_matrix(index,:) = [];
magnitude(index) = [];
Rjb(index,:) = [];
Vs(index,:) = [];
mechanism(index,:) = [];


SS = zeros(length(mechanism),1);
N = zeros(length(mechanism),1);
R = zeros(length(mechanism),1);
U = zeros(length(mechanism),1);

SS(mechanism == 0) = 1;
N(mechanism == 1) = 1;
R(mechanism == 2) = 1;
U(mechanism > 2) = 1;

mechanism_BA = zeros(length(mechanism),1);
mechanism_BA(mechanism > 2) = 1;
mechanism_BA(mechanism == 0) = 2;
mechanism_BA(mechanism == 1) = 3;
mechanism_BA(mechanism == 2) = 4;

end