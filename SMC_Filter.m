%**************************************************************************
%     Function: filtro
%     Last modified: 29/09/16
%     Author: Juan D. Cardenas
%     Description: this algorithm computes a filtered state via artificial
%                  state generation, using the particle filter method
%     Input:
%           * X_filter_K_1: Last filtered system state
%           * meas: Current real system measurement
%           * u: Current system input
%           * step: step time [integration purposes]
%           * Np: Number of particles
%           * e_k_p: Current e value [ARIMA noise]
%
%     Output:
%           * X_filter: Filtered state
%           * wk: Current weight vector
%           * e_k_p: Current e value [ARIMA noise]
%**************************************************************************
function [X_filter, wk, e_k_p] = SMC_Filter(x_filter_K_1, meas, u, wk, step, Np, e_k_p)

% Memory allocation
	x_k = zeros(6, Np);                                                     % Particles memory allocation
    pdf_y_x = zeros(1, Np);                                                 % Particles memory allocation

% Get past weights
    w_k_1 = wk;                                                             % Weight assignation

% System noise parameters with Gaussian assumption
    mu_w = 0;                                                               % System noise mean for particles
    sigma_w = sqrt(10);                                                     % System noise variance for particles
    % Noise generation
    noise_w = normrnd(mu_w,sigma_w, 7, Np);                                 % System noise generation
    [noise_v, e_k_p] = noise_ARIMA(1, e_k_p, Np);                           % Measurement noise generation with ARIMA assumption
    
    for i = 1:Np                                                            % Artificial state (particles) generation
        
        dot_particles_k = modelNL_filtro(x_filter_K_1, u);                  % Compute current state derivates
        X_k_aux = dot_particles_k.*step + x_filter_K_1;                     % Get current state via integration (euler method)
        x_k([1 3:6], i) = X_k_aux([1 3:6]) + noise_w([1 3:6], i);           % System noise addition
        
        x_k(2, i) = X_k_aux(2) + noise_v(i);                                   % Measurement noise addition
        
        pdf_y_x(i) = pdf_fun(meas, x_k(:, i));                              % Computes the pdf(y|x)
    end
    
    wk = bsxfun(@times, w_k_1, pdf_y_x);                                    % Computes current weights
	wk = bsxfun(@rdivide, wk, sum(wk));                                     % Normalize current weights
    
	Neff = 1/sum(wk.^2, 2);                                                 % Computes how many particles are efective
	Nt = 0.5*Np;                                                            % 50% of the particles must be effective otherwise, they should be resampled
    if Neff < Nt                                                            % Resampling function
        [x_k, wk] = resampling(x_k, wk);
    end

    pre_state = bsxfun(@times, x_k, wk);                                    % Computes the importance of each particle
    X_filter = sum(pre_state, 2);                                           % Computes the filtered state
end

%**************************************************************************
%     Function: pdf_fun
%     Last modified: 29/09/16
%     Author: Juan D. Cardenas
%     Description: this algorithm computes the probability function p(y_k|x_k)
%     Input:
%           * Z: Real measurement value
%           * Z_hat: Artificial measurement value
%
%     Output:
%           * p_y_x: p(y_k|x_k) value for the Artificial measurement value
%**************************************************************************
function p_y_x = pdf_fun(Z, Z_hat)
    mu_v = 0;                                                               % Measurement for real noise mean
    sigma_v = sqrt(1);                                                    % Measurement for real noise variance
    
	z_m = Z_hat(2);                                                         % Artificial measurement
    
    error = (Z - z_m);                                                      % Computes an error between real and artificial measurements
    
    lamnda = 1/sqrt(det(sigma_v)*(2*pi)^(1));                               % Auxiliar value for p(y|x)
    p_y_x = lamnda*exp(-0.5*(error'/sigma_v)*error);                        % Probability function p(y|x)
end

%**************************************************************************
%     Function: resampling
%     Last modified: 29/09/16
%     Author: Juan D. Cardenas
%     Description: this algorithm resamples the particles in order to avoid
%     degenation process in the particle filter
%     Input:
%           * xk: Particles
%           * wk: Particles' weights
%
%     Output:
%           * x_h: Resampled particles
%           * w_h: Resampled weights
%**************************************************************************
function [x_h, w_h] = resampling(xk, wk)

    Np = length(wk);                                                        % Get the number of particles
    W = cumsum(wk);                                                         % Generates a cumulative distribution function (CDF) using the weights 

    index = zeros(1, Np);                                                   % Memory allocation, vector with selected particles' indexes
    j = 1; i = 1;                                                           % Index inicialization
    
    while j <= Np                                                           % Roulette algorithm 
       ruleta = rand();                                                     % Generates a uniform random number
       if ruleta <= W(i)                                                    % If the CDF evaluated in the i^th particle is bigger than the random number, the particle is selected in the resampling process                                                
            index(j) = i;                                                   % Save particle's index
            j = j + 1;                                                      % Next iteration
       end                                                                  % The algorithm ends when there are Np number of selected particles
       
       i =  round(Np*rand() + 1);                                           % Generates an aleatory particle's index
       
       if i > Np                                                            % Only index smaller or equal to Np are accepted 
           i = Np;
       end
    end
    
    x_h = xk(:, index);                                                     % Save selected particles in the resampling process
    w_h_aux = W(index);
	w_h = bsxfun(@rdivide, w_h_aux, sum(w_h_aux));
end