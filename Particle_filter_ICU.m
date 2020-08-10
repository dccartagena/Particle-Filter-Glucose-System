%**************************************************************************
%     Particle filter implementation in a Insulin-Nutrition-Glucose 
%     (ICING) Model (Lin, J. et al (2011)) 
%     Last modified: 29/09/16
%     Author: Jose Fernando Garcia, Estefania Aguirre Zapata 
%     and Juan D. Cardenas
%**************************************************************************

%% Inicialization
close all;
clear all;
clc;

%% Run simulink simulation
ModelParam;
model = 'Model_7'; 
load_system(model);
sim(model);

%% Auxiliar functions
% Noise with gaussian assumption
mu_w = 0; % State noise mean
sigma_w = 10; % State noise variance 
noise_w = @(m) normrnd(mu_w,sigma_w, m, 6);                                 % State noise function

%% Simulation parameter
step = tiempo_sim(2) - tiempo_sim(1);                                       % Time step
iter = length(tiempo_sim);                                                  % Number of iteration

%% Memory allocation
v = zeros(iter, 1);                                                         % Measurement noise

%% Parameter inicialization
% X(1): Blood glucose [mmol/L]
% X(2): Interstitial glucose [mmol/L]
% X(3): Interstitial insulin [mU/L]
% X(4): Plasmatic insulin [mU/L]
% X(5): Glucose in the stomach [mmol]
% X(6): Glucose in the gut [mmol]

X_0 = [5 5 74 108.58 0 0]';                                                 % First state

%% Noise generation                                                         
mu = zeros(1, 6);
sigma = [1 0 0      0      0  0;
         0 1 0      0      0  0;
         0 0 0.0165 0      0  0;
         0 0 0      0.0142 0  0;
         0 0 0      0      11 0;
         0 0 0      0      0  16];
w = mvnrnd(mu, sigma, iter);                                                % System noise generation
     
% Noise generation with ARIMA assumption
e_k_1 = 0;
for i = 1:iter
    [v(i), e_k_1] = noise_ARIMA(i, e_k_1, 1);
end

X_noise(:, [1 3:6]) = Xreal(:, [1 3:6]) + w(:, [1 3:6]);                    % Noise system addition
X_noise(:, 2) = Xreal(:, 2) + v;                                            % Measurements and noise addition
meas = X_noise(:, 2);

%% Inputs
% u_1: Insulin [mU/min]
% u_2: Intravenous feed [mU/min]
% u_3: Enteral feed [mmol/min]

u = [u_1 u_2 u_3];                                                          % Inputs from simulink simulation
%% Particle filter parameters
Np = 5000;                                                                   % Number of particles

X_filter = zeros(6, iter);                                                  % Memory allocation for filtered states

X_filter(:, 1) = X_noise(1, :)';                                            % First state with system noise

%% Particle filter implementation
tic;
wk = repmat(1/Np, [1, Np]);                                                 % Weight initialization
e_k_p = 0;
for i = 2:iter   
    [X_filter(:, i), wk, e_k_p] = SMC_Filter(X_filter(:, i-1), meas(i), ...     % Filtered state
                                         u(i, :), wk, step, Np, e_k_p);    
end
tiempo = toc

delta_X = Xreal' - X_filter;                                                % Error computation
RMSE_aux = sqrt(sum(delta_X .^ 2, 2) / iter);                        

RMSE = sum(RMSE_aux)

%% Figures and plots comparing: real state, filtered state, ideal state

% X(1): Blood glucose [mmol/L]
figure; 
plot(tiempo_sim, X_noise(:, 1), 'b', tiempo_sim, X_filter(1, :), 'r', tiempo_sim, Xreal(:, 1), 'g');
title('Blood glucose'); 
legend('Real state','Filtered state', 'Ideal state');
xlabel('time [s]'); ylabel('Glucose concentration [mmol/l]');
grid on; grid minor;

% X(2): Interstitial glucose [mmol/L]
figure; 
plot(tiempo_sim, X_noise(:, 2), 'b', tiempo_sim, X_filter(2, :), 'r', tiempo_sim, Xreal(:, 2), 'g');
title('Interstitial glucose - measurement')
legend('Real state','Filtered state', 'Ideal state');
xlabel('time [s]'); ylabel('Glucose concentration [mmol/l]');
grid on; grid minor;

% X(3): Interstitial insulin [mU/L]
figure; 
plot(tiempo_sim, X_noise(:, 3), 'b', tiempo_sim, X_filter(3, :), 'r', tiempo_sim, Xreal(:, 3), 'g');
title('Interstitial insulin'); 
legend('Real state','Filtered state', 'Ideal state');
xlabel('time [s]'); ylabel('Insulin concentration [mU/l]');
grid on; grid minor;

% X(4): Plasmatic insulin [mU/L]
figure; 
plot(tiempo_sim, X_noise(:, 4), 'b', tiempo_sim, X_filter(4, :), 'r', tiempo_sim, Xreal(:, 4), 'g');
title('Plasmatic insulin'); 
legend('Real state','Filtered state', 'Ideal state');
xlabel('time [s]'); ylabel('Insulin concentration [mU/l]');
grid on; grid minor;

% X(5): Glucose in the stomach [mmol]
figure; 
plot(tiempo_sim, X_noise(:, 5), 'b', tiempo_sim, X_filter(5, :), 'r', tiempo_sim, Xreal(:, 5), 'g');
title('Glucose in the stomach'); 
legend('Real state','Filtered state', 'Ideal state');
xlabel('time [s]'); ylabel('Glucose [mmol]');
grid on; grid minor;

% X(6): Glucose in the gut [mmol]
figure; 
plot(tiempo_sim, X_noise(:, 6), 'b', tiempo_sim, X_filter(6, :), 'r', tiempo_sim, Xreal(:, 6), 'g');
title('Glucose in the gut'); 
legend('Real state','Filtered state', 'Ideal state');
xlabel('time [s]'); ylabel('Glucose [mmol]');
grid on; grid minor;
