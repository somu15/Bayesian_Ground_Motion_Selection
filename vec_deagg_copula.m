%% Performs vector deaggregation given scalar deaggregations
%  This code cannot be used for more than two IMs

% "A Bayesian Treatment of the Conditional Spectrum Approach for Ground
% Motion Selection". Report by Somayajulu Dhulipala, Jack Baker and
% Madeleine Flint.
%
% Created by Somayajulu Dhulipala, 20 May 2017
% Modified by Somayajulu Dhulipala, 23 December 2017

%% INPUTS

% Deagg            =          matrix consisting of scalar deaggregations
% (NM X NR X N_IMs)

% NM               =          number of magnitude bins
% NR               =          number of distance bins
% N_IMs            =          2;
% Deagg_min        =          deaggregation matrix for a low value of IM
% (NM X NR)

% Haz_matrix       =          hazard matrix consisting of required hazard
% values in the first column and the hazard value for a low IM level in the
% second column (N_IMs X 2)

% corr_req         =          Pearson correaltion between the two IMs 
% plots_req        =          'Yes' or 'No'

%% OUTPUTS

% deagg_final      =          vector deaggregation

%% EXPLICITLY DEFINED FUNCTIONS/ROUTINES USED

% NONE

%% FUNCTION BEGINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [deagg_final] = vec_deagg_copula(Deagg, Deagg_min, Haz_matrix, corr_req, N_IMs, NM, NR, plots_req)

for ii = 1:NM
    for jj = 1:NR
        Z_vec = zeros(N_IMs,1);
        for kk = 1:N_IMs
            temp = Deagg(ii,jj,kk)*Haz_matrix(kk,1);
            temp_min(kk) = Deagg_min(ii,jj)*Haz_matrix(kk,2);
            CDF_norm = 1-temp/temp_min(kk);
            CDF_norm(isnan(CDF_norm)) = 0;
            CDF_norm(isinf(CDF_norm)) = 0;
            sto_CDF(kk) = CDF_norm;
            Z_vec(kk) = norminv(CDF_norm);
        end
        Z_vec = reshape(Z_vec, N_IMs,1);
        CDF_unif = mvncdf(Z_vec,[0;0],corr_req);
        
        CCDF = 1 - (sto_CDF(1)) - (sto_CDF(2)) + CDF_unif;
        
        CCDF(isnan(CCDF)) = 0;
        CCDF(isinf(CCDF)) = 0;
        CCDF = CCDF*prod(temp_min);
        deagg_final(ii,jj) = CCDF;
    end
end

deagg_final = deagg_final/sum(sum(deagg_final));

if strcmp(plots_req, 'Yes') == 1
    bar3(deagg_final)
end
end