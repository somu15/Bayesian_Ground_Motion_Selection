%% Prepares the input file for EXSIM accelerogram simulations
%
% "A Bayesian Treatment of the Conditional Spectrum Approach for Ground
% Motion Selection". Report by Somayajulu Dhulipala, Jack Baker and
% Madeleine Flint.
%
% Created by Somayajulu Dhulipala, 28 August 2017
% Modified by Somayajulu Dhulipala, 22 December 2017

%% INPUTS

file_name    =       'C:\...'; % the current directory 
% (the above directory should contain the EXSIM program and 
% EXSIM_simulations_inputs.m routine)
Niter        =            500; % number of ground motion simulations

%% OUTPUTS

% param(:,1) = Mw; param(:,2) = stress drop in bars; param(:,3) = cood1 in
% rads; param(:,4) = cood2 in rads; param(:,5) = Fault Type; param(:,6) =
% strike; param(:,7) = dip;

% per        =        Time Period (s)
% Sa         =        Spectral Acc. (g)

%% EXPLICITLY DEFINED FUNCTIONS/ROUTINES USED

% 1) EXSIM_simulations_inputs.m
% 2) EXSIM12.exe (EXSIM program)

%% USER INVOLVEMENT ENDS HERE
%% INITIALIZE VARIABLES

cd(filename);

%% PERFORM EXSIM SIMULATIONS AND EXTRACT OUTPUTS

for kkk = 1:Niter
    
EXSIM_simulations_inputs

!EXSIM12.exe Inputs.txt

max_ind = 1;
base_name1 = 'PSA\M6d5S130S001iter00';

param(kkk,1) = Mw; param(kkk,2) = SD; param(kkk,3) = cood1; 
param(kkk,4) = cood2; param(kkk,5) = Fault; param(kkk,6) = strk;
param(kkk,7) = dip;

fid = fopen(strcat(base_name1,num2str(1),'.tmp'));
    
for ii = 1:14
    temp = fgetl(fid);
    if ii == 10
        param(kkk,8) = str2double(temp(20:length(temp)));
    end
end

temp = str2num(fgetl(fid));
sto_ind = 0;

while temp ~= -1
    sto_ind = sto_ind+1;
    if kkk == 1
    per(sto_ind) = temp(2);
    end
    Sa(sto_ind,kkk) = temp(4)*0.00101;
    temp = str2num(fgetl(fid));
end

fclose all
progressbar(kkk/Niter)

end

clearvars -except param per Sa

figure
loglog(per, Sa)
grid on
xlabel('Time Period (s)')
ylabel('Sa (g)')
