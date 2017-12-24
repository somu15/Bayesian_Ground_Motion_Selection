%% Extracts the prior distributions for beta, Delta and Q

% "A Bayesian Treatment of the Conditional Spectrum Approach for Ground
% Motion Selection". Report by Somayajulu Dhulipala, Jack Baker and
% Madeleine Flint.
%
% Created by Somayajulu Dhulipala, 20 May 2017
% Modified by Somayajulu Dhulipala, 23 December 2017

% PLEASE REPORT ANY ISSUES OR PROBLEMS WITH THIS CODE TO
% lakshd5@vt.edu

%% OUTPUTS

% CMS          =          the Conditional Mean Spectrum
% TV           =          the Target Variability

%% EXPLICITLY DEFINED FUNCTIONS/ROUTINES USED

% 1) BA_X.m
% 2) BA_X_Large.m
% 3) CleanData_V7.m
% 4) pinky.m
% 5) Prior_Extraction.m
% 6) vec_deagg_copula.m
% 7) progressbar.m
% 8) Number_Finder_V1.m

%% You can directly run this code by clicking the play button

%% User inputs

clearvars;
clc;

prompt = 'Do you have a ready to use, Bayesian amenable workspace? (0-No, 1-Yes)';
decision1 = input(prompt);

if decision1 == 0

prompt = 'Enter the current directory in a string';
file_name = input(prompt);
cd(file_name)

%Oscillator periods

display('These are oscillator periods considered for the analysis:')
req_per = [0.100000000000000,0.120000000000000,0.140000000000000,0.160000000000000,0.190000000000000,0.220000000000000,0.260000000000000,0.300000000000000,0.360000000000000,0.420000000000000,0.480000000000000,0.550000000000000,0.667000000000000,0.800000000000000,0.900000000000000,1.10000000000000,1.30000000000000,1.50000000000000,1.70000000000000,2,2.40000000000000,2.80000000000000,3.20000000000000,3.80000000000000,4.60000000000000,5.50000000000000,6,7.50000000000000,8.50000000000000,10];
disp(req_per')
prompt = 'Do you want to make changes? (1 = Yes; 0 = No)';
decision = input(prompt);
if decision == 1
    prompt = 'Enter the oscillator period vector';
    req_per = input(prompt);
end

% Ground Motion Database

prompt = 'Do you want to use the default NGA West 2 database for analysis? (1 = Yes; 0 = No)';
decision = input(prompt);
if decision == 0
    prompt = 'Enter the file path for spectral acceleration matrix. It should be a .mat file with N_o observations at N_t required oscillator periods.';
    file = input(prompt);
    K = importdata(file);
    prompt = 'Enter the file path for magnitude vector. It should be a .mat file with N_o observations.';
    file = input(prompt);
    magnitude = importdata(file);
    prompt = 'Enter the file path for Joyner-Boore distance (Km.) vector. It should be a .mat file with N_o observations.';
    file = input(prompt);
    Rjb = importdata(file);
    prompt = 'Enter the file path for Soil Vs_30 (m/s) vector. It should be a .mat file with N_o observations.';
    file = input(prompt);
    Vs = importdata(file);
    prompt = 'Enter the file path for fault mechanism vector. It should be a .mat file with N_o observations. The following is a description of mechanisms: 0 - Strike slip; 1 - Normal; 2 - Reverse; 3 - Others.';
    file = input(prompt);
    mechanism = importdata(file);
    prompt = 'Enter the number of non-spectral IMs to consider.';
    N_IMs = input(prompt);
    for ii = 1:N_IMs
        prompt = strcat('Index of the non-spectral IM is-',num2str(ii),'.',' Enter IM',num2str(ii),' vector. It should be N_o observations.');
        temp = input(prompt);
        K = [K temp];
    end
    IM_matrix = K;
else
    display('Two non-spectral IMs are considered: PGA (index 1), PGV (index 2).')
    N_IMs = 2;
    load('Final_data_cleaned.mat');
    dummy = 1;
    for ii = 1:(length(req_per)+N_IMs)
    if ii <= length(req_per)
index(ii) = find(Periods == req_per(ii));
    else
        index(ii) = length(Periods)+dummy;
        dummy = dummy + 1;
    end
    end
    IM_matrix = IM_matrix(:,index);
end

[IM_matrix, magnitude, Rjb, Vs, mechanism,SS,N,R,U] = CleanData_V7(IM_matrix, magnitude, Rjb, Vs, mechanism, req_per, [0.01 2;1 250]);
M_ref = 4.5; R_ref = 1; M_h = 6.5; Vs_ref = 760; h = 2.1668;

[X] = BA_X(magnitude, Rjb, Vs, SS, N, R, U, M_ref, R_ref, M_h, Vs_ref, h);

Y = log(IM_matrix);

clearvars -except Y X req_per magnitude M_ref R_ref M_h Vs_ref h N_IMs

[N_o, N_t, N_p, N_e, N_re] = Number_Finder_V1(Y, magnitude,min(size(X)), N_IMs);

% Conditioning IMs

prompt = 'Enter the number of conditioning IMs';
N_cond = input(prompt);

prompt = 'Of which the number of non-spectral IMs are';
N_Ncond = input(prompt);

if (N_cond-N_Ncond) > 0
        prompt = 'Enter the vector of oscillator periods of conditioned spectral IMs';
        per_deagg = input(prompt);
        ind_vec = find(ismember(req_per, per_deagg));
else
    ind_vec = [];
end



if N_Ncond > 0
    for ii = 1:N_Ncond
        prompt = 'Enter the index of non-spectral IM';
        temp = input(prompt);
        ind_vec = [ind_vec (temp+N_t)];
    end
end


prompt = 'Input the IM level vector corresponding to the hazard levels for all the conditioning IMs (spectral IMs should come first and then non spectral IMs according to their indices)';
IM_vec = input(prompt);

ss = 0;
if (N_cond-N_Ncond) > 0
        for ii = 1:(N_cond-N_Ncond)
            ss = ss + 1;
            prompt = 'Input the hazard values for spectral IM at design level and minimum level [design_level_hazard    minimum_IM_level_hazard]';
            Haz_matrix(ss,:) = input(prompt);
            prompt = 'Input Magnitude and Distance deaggregation information for spectral IM at the corresponding hazard level [a deaggregation matrix]';
            temp = input(prompt);
            Deagg(:,:,ss) = temp;
            
        end
end
            
if N_Ncond > 0
        for ii = 1:N_Ncond
            ss = ss + 1;
            prompt = 'Input the hazard values for non-spectral IM at design level and minimum level [design_level_hazard    minimum_IM_level_hazard]';
            Haz_matrix(ss,:) = input(prompt);
            prompt = 'Input Magnitude and Distance deaggregation information for non-spectral IM at the corresponding hazard level [a deaggregation matrix]';
            temp = input(prompt);
            Deagg(:,:,ss) = temp;
        end
end

prompt = 'Input Magnitude and Distance deaggregation information for spectral IM at its minimum level [a deaggregation matrix]';
            temp = input(prompt);
            Deagg_min(:,:) = temp;
            
prompt = 'Enter the Magnitude bins vector';
M = input(prompt);
NM = length(M);

prompt = 'Enter the Distance bins vector';
R = input(prompt);
NR = length(R);

prompt = 'Input the site surface shear wave velocity (Vs30 m/s)';
VS = input(prompt);

% Analysis inputs

prompt = 'Enter the total number of iterations';
N_iter = input(prompt);
prompt = 'Enter the iterations until burn-in';
N_burn = input(prompt);
prompt = 'Enter the iterations for vector deaggregation';
N_deagg = input(prompt);

prompt = 'Do you want to use EXSIM priors? (1-YES, 0-NO)';
prior_ID = input(prompt);

else
    prompt = 'Enter the .mat file name containing Bayesian amenable workspace';
    mat_file = input(prompt);
    load(mat_file)
    prompt = 'Do you want to use EXSIM priors? (1-YES, 0-NO)';
    prior_ID = input(prompt); 
end

%% Main code

IM_vec = reshape(IM_vec,N_cond,1);
deagg_count = 1;
temp = size(Y);
N_t_final = temp(2);

alpha = 5*ones(N_t_final, N_p);
a = 2; b = 5;
PHI = eye(max(size(Y)));
sss = 1;
SIGMA = eye(N_t_final);
pp = 1;

if prior_ID == 1
[beta_0,DELTA,Q] = Prior_Extraction(ID);
else
beta_0 = zeros(N_t_final*N_p,1);
DELTA = 10*eye(N_t_final*N_p);
Q = 10*eye(N_t_final);
end

nu = 31+(N_t_final-length(req_per));

beta_cap = reshape(Y'*X*inv(X'*X),N_t_final*N_p,1);



for ii = 1:N_iter
    
    Y_beta = normrnd(0,1,N_t_final*N_p,1);
    Y_sigma = normrnd(0,1,N_o+N_t_final+nu+1,N_t_final);
    
    P1 = kron(X'*X, inv(SIGMA));
    A_beta = chol(inv(inv(DELTA)+P1))';
    
    P2 = Y-X*alpha';
    A_sigma = chol(P2'*P2+Q)';

M_beta = inv(inv(DELTA)+P1)*(inv(DELTA)*beta_0+P1*beta_cap);
    
beta = reshape(M_beta, N_t_final, N_p);

alpha = A_beta*Y_beta+M_beta;
alpha = reshape(alpha, N_t_final, N_p);
SIGMA = A_sigma*inv(Y_sigma'*Y_sigma)*A_sigma';

sto_alpha(:,:,ii) = alpha;
sto_SIGMA(:,:,ii) = SIGMA;

if ii>N_burn && ii<=(N_burn+N_deagg)
    sto_alpha(:,:,deagg_count) = alpha;
    sto_SIGMA(:,:,deagg_count) = SIGMA;
    deagg_count = deagg_count+1;
    
    if ii==(N_burn+N_deagg)
        for kk = 1:N_t_final
    for jj = 1:N_p
        deagg_mean_alpha(kk,jj) = mean(sto_alpha(kk,jj,:));
    end
        end

for kk = 1:N_t_final
    for jj = 1:N_t_final
        deagg_mean_SIGMA(kk,jj) = mean(sto_SIGMA(kk,jj,:));
    end
end

if N_cond > 1
    [corr_req] = compute_req_corr(deagg_mean_SIGMA, ind_vec);
    fMR12 = vec_deagg_copula(Deagg, Deagg_min, Haz_matrix, corr_req, N_cond, NM, NR, 'No');
else
    fMR12 = Deagg;
end

for ss = 1:NM
    temp_M(ss) = sum(fMR12(ss,:));
end
temp_M = temp_M/sum(temp_M);
M_PSHA = sum(M.*temp_M');
for ss = 1:NR
    temp_R(ss) = sum(fMR12(:,ss));
end
temp_R = temp_R/sum(temp_R);
R_PSHA = sum(R.*temp_R');

    end
end

    if ii>(N_burn+N_deagg) 
        
       M_sto(pp) = M_PSHA;
     R_sto(pp) = R_PSHA;
     pp = pp + 1;
       fault = randi(3)-1;
%        
if ID > 0
    MU = alpha*BA_X(M_PSHA, R_PSHA, VS,  0, 0, 0, 1, M_ref, R_ref, M_h, Vs_ref, h)';
else
MU = alpha*BA_X_Large_M(M_PSHA, R_PSHA, VS,  0, 0, 0, 1, M_ref, R_ref, M_h, Vs_ref, h)';
end
        MU_2 = MU(ind_vec);
       MU_1 = MU;
       MU_1(ind_vec) = [];
       
       [SIG_11, SIG_22, SIG_12] = Mat_Trans_V1(SIGMA, ind_vec, N_t_final);
       
       MU_COND = MU_1 + SIG_12*inv(SIG_22)*(log(IM_vec)-MU_2);
       SIG_COND = SIG_11 - SIG_12*inv(SIG_22)*SIG_12';
       
       SIG_COND = triu(SIG_COND)+triu(SIG_COND)'-diag(diag(SIG_COND));
       
       Y_pred = (mvnrnd(MU_COND,SIG_COND));
       
       Y_pred = rearr_Y(Y_pred, ind_vec, IM_vec);
       
       sto(sss,:) = Y_pred;
       
       sss = sss+1;

    end

    progressbar(ii/N_iter)
    
end

K = sto_alpha(:,:,1001:N_iter);
K1 = sto_SIGMA(:,:,1001:N_iter);
for ii = 1:N_t_final
    for jj = 1:N_p
        mean_alpha(ii,jj) = mean(K(ii,jj,:));
        std_alpha(ii,jj) = std(K(ii,jj,:));
    end
end

for ii = 1:N_t_final
    for jj = 1:N_t_final
        mean_SIGMA(ii,jj) = mean(K1(ii,jj,:));
        std_SIGMA(ii,jj) = std(K1(ii,jj,:));
    end
end

Y_pred = X*mean_alpha';

std_req = sqrt(diag((mean_SIGMA))); 
if ID > 0
    MU_mine = mean_alpha*BA_X(M_PSHA, R_PSHA, VS,  0, 0, 0, 1, M_ref, R_ref, M_h, Vs_ref, h)';
else
MU_mine = mean_alpha*BA_X_Large_M(M_PSHA, R_PSHA, VS,  0, 0, 0, 1, M_ref, R_ref, M_h, Vs_ref, h)';
end

eps_req = (log(IM_vec)-MU_mine(ind_vec))/std_req(ind_vec);


CMS = exp(mean(sto));
TV = std((sto));

clearvars -except req_per CMS TV
 
figure
loglog(req_per, (CMS(1:N_t)),'linewidth',1.5);
xlabel('Time Period (s)')
ylabel('SA (g)')
title('CMS')
% 
figure
semilogx(req_per, TV(1:N_t), 'linewidth',1.5)
xlabel('Time Period (s)')
ylabel('Conditional STD')
title('TV')
