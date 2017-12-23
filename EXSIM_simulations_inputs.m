%% Prepares the input file for EXSIM accelerogram simulations
%
% "A Bayesian Treatment of the Conditional Spectrum Approach for Ground
% Motion Selection". Report by Somayajulu Dhulipala, Jack Baker and
% Madeleine Flint.
%
% Created by Somayajulu Dhulipala, 28 August 2017
% Modified by Somayajulu Dhulipala, 22 December 2017

%% INPUTS

file_name    =  'C:\...'; % the current directory 
% (the above directory should contain an EXSIM "Input.txt" file)
Mw_high      =       9.5; % largest moment magnitude
Mw_low       =       7.0; % smallest moment magnitude 
SD_low       =      50.0; % smallest static stress drop in bars
SD_high      =     150.0; % largest static stress drop in bars
cood1_low    =      0.01; % smallest coordinate 1 from origin in rads
cood1_high   =      0.31; % largest coordinate 1 from origin in rads
cood2_low    =      0.01; % smallest coordinate 2 from origin in rads
cood2_high   =      0.31; % largest coordinate 2 from origin in rads
strk         =       0.0; % strike of the fault in degrees
dip          =      90.0; % dip of the fault in degrees
Fault        =  randi(4); % FT randomized b/w SS, N, R, and U

%% OUTPUTS

% There are no actual outputs associated with this code rather than setting
% up the input file for EXSIM ground motion simulations in a folder of the
% user's choice.

%% EXPLICITLY DEFINED MATLAB FUNCTIONS/ROUTINES USED

% NONE

%% USER INVOLVEMENT ENDS HERE
%% INITIALIZE VARIABLES

cd(file_name)
Mw           =        Mw_low + (Mw_high - Mw_low)*rand;
SD           =        SD_low + (SD_high - SD_low)*rand;
cood1        =        cood1_low + (cood1_high - cood1_low)*rand; 
cood2        =        cood2_low + (cood2_high - cood2_low)*rand;

%% PREPARE THE INPUT EXSIM FILE

fid = fopen('Inputs.txt','r');
i = 1;
tline = fgetl(fid);
A{i} = tline;

while ischar(tline)
    i = i+1;
    tline = (fgetl(fid));
    A{i} = tline;
end
fclose(fid);

A = A';
A{8} = char(strcat({' '},num2str(Mw),{' '},num2str(SD),{' 1 0.035'}));
A{12} = char(strcat({' '},num2str(strk),{' '},num2str(dip),{' 3.0'}));
A{114} = char(strcat({' '},num2str(cood1),{' '},num2str(cood2)));

if Mw>8
    A{30} = ' 350.0 0.0 2.0 2.0 70.0 !fault length and width, dl, dw, stress_ref';
else
    A{30} = ' 0.0 0.0 2.0 2.0 70.0 !fault length and width, dl, dw, stress_ref';
end

if Fault == 1
    A{15} = ' S';
elseif Fault == 2
    A{15} = ' R';
elseif Fault == 3
    A{15} = ' N';
else
    A{15} = ' U';
end

fid = fopen('Inputs_Test.txt', 'w');
for i = 1:numel(A)
    if strcmp(A{i},-1) == 1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end

fclose(fid);
