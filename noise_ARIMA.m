%**************************************************************************
%     Function: noise_ARIMA
%     Last modified: 29/09/16
%     Author: Juan D. Cardenas
%     Description: this algorithm creates artificial noise from an ARMA
%     method
%     Input:
%           * t: Current iteration
%           * e_k_1: Last e value
%
%     Output:
%           * Error: Noise generated via ARIMA method
%           * e_k: Current e value
%**************************************************************************
function [Error, e_k] = noise_ARIMA(t, e_k_1, n)
    % Parameters
    lamnda = 15.96;
    epsilon = -5.471;
    delta = 1.6898;
    gamma = -0.5444;
    
    % Gaussian noise parameters
    mu = 0;                                                                 % Mean of random number for noise generation
    sigma = sqrt(1);                                                      % Standard deviation of random number for noise generation
    v_k = randn(1, n).*sigma + mu;                                              % Random number generation with gaussian assumption                          
    
    if t == 0                                                               % The e value depends on the iteration
        % ARIMA noise with assumption
        e_k = v_k;                                                          % e value in the first iteration
    else
        e_k = 0.7.*(e_k_1 + v_k);                                           % e value in the next iterations
    end
    
    Error = epsilon + lamnda.*sinh((e_k - gamma)./(delta));                 % Noise computation
end