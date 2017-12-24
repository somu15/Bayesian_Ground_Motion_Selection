%% Extracts the prior distributions for beta, Delta and Q

% "A Bayesian Treatment of the Conditional Spectrum Approach for Ground
% Motion Selection". Report by Somayajulu Dhulipala, Jack Baker and
% Madeleine Flint.
%
% Created by Somayajulu Dhulipala, 20 May 2017
% Modified by Somayajulu Dhulipala, 23 December 2017

%% INPUTS

% ID        =         0 when all the records have a Mw > 7
%           =         1 otherwise

% Simulated records database in a manner amenable to Bayesian analysis. To
% achieve this, run the first part (named as ???) in the Master.m code,
% follow the instructions and then save the output variables as For_Priors.mat

%% OUTPUTS

% beta_0    =         the prior mean coefficient matrix ((Nt+Nn) X 11(9))
% Nt is the number of time periods, Nn is the number of non-spectral IMs
% 11 indicates both small and large Mw earthquakes have been considered
% 9 indicates only large Mw earthquakes have been considered

% DELTA     =         the covariance matrix on the GMPM coefficients
% (11(9) X 11(9))

% Q         =         Scale matrix for the covariance matrix on GMPM
% residuals ((Nt+Nn) X (Nt+Nn))

%% EXPLICITLY DEFINED FUNCTIONS/ROUTINES USED

% NONE

%% FUNCTION BEGINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [beta_0, DELTA, Q] = Prior_Extraction(ID)

load('EXSIM_Priors1.mat')

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
beta_0 = zeros(N_t_final*N_p,1);
DELTA = 10*eye(N_t_final*N_p);
nu = 31+(N_t_final-length(req_per));
Q = 10*eye(N_t_final);
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

if ID > 0
    
sto_alpha(:,:,ii) = [alpha(:,1:4) normrnd(0,sqrt(10),30,1)...
normrnd(0,sqrt(10),30,1) alpha(:,5:9)];

else 
    
sto_alpha(:,:,ii) = alpha;

end

sto_SIGMA(:,:,ii) = SIGMA;

progressbar(ii/N_iter)
    
end

K = sto_alpha(:,:,1001:N_iter);
K1 = sto_SIGMA(:,:,1001:N_iter);

if ID > 0
    
K = reshape(K,N_t_final*(N_p+2),2000);

else
    
    K = reshape(K,N_t_final*(N_p),2000);
    
end

mean_alpha = mean(K')';
cov_alpha = cov(K');

for ii = 1:N_t_final
    for jj = 1:N_p
        
        mean_alpha1(ii,jj) = mean(K(ii,jj,:));
        var_alpha(ii,jj) = var(K(ii,jj,:));
        
    end
end

for ii = 1:N_t_final
    for jj = 1:N_t_final
        
        scale_SIGMA(ii,jj) = var(K1(ii,jj,:));
        
    end
end

beta_0 = mean_alpha; DELTA = cov_alpha; Q = scale_SIGMA;
 