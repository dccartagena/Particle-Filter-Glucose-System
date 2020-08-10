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
sigma_w = 0;%0.1; % State noise variance 
noise_w = @(m) normrnd(mu_w,sigma_w, m, 6);                                 % State noise function

mu_v = 0;                                                                   % Measurement noise mean
sigma_v = sqrt(5);                                                          % Measurement noise variance
noise_v = @(m) normrnd(mu_v,sigma_v, m, 1);                                 % Measurement noise function

% ARIMA noise with assumption
mu_ARIMA = 0;                                                               % Measurement noise mean with ARIMA assumption
sigma_ARIMA = sqrt(0.1);                                                    % Measurement noise variance with ARIMA assumption
e_k_1 = randn().*sigma_ARIMA + mu_ARIMA;                                    % First e value in ARIMA assumption

%% Simulation parameter
iter = length(tiempo_sim);                                                  % Number of iteration for filter simulation

%% Memory allocation
v = zeros(iter, 1);                                                         % Measurement noise
meas = zeros(1, iter);                                                      % Real measurement

%% Parameter inicialization
% X(1): Blood glucose [mmol/L]
% X(2): Interstitial glucose [mmol/L]
% X(3): Interstitial insulin [mU/L]
% X(4): Plasmatic insulin [mU/L]
% X(5): Glucose in the stomach [mmol]
% X(6): Glucose in the gut [mmol]

X_0 = [5.5 5.5 81.40 119.44 0 0 0]';                                        %First state case 1
%X_0 = [4.5 4.5 66.6 97.72 0 0 0];                                          %First state case 2
%X_0 = [6 6 88.8 130.30 0 0 0];                                               %First state case 3
%X_0 = [4 4 59.2 86.86 0 0 0];                                                %First state case 4
%X_0 = [6.5 6.5 96.2 141.15 0 0 0];                                           %First state case 5
%X_0 = [3.5 3.5 51.8 76 0 0 0];                                               %First state case 6

%% Noise generation
w = noise_w(iter);                                                          % System noise generation
%v(:) = noise_v(iter);                                                      % Measurement noise generation with Gaussian assumption

for i = 1:iter
    [v(i), e_k_1] = noise_ARIMA(i, e_k_1);                                  % Noise generation with ARIMA assumption
end

X_noise(:, [1 3:6]) = Xreal(:, [1 3:6]) + w(:, [1 3:6]);                    % Noise system addition

meas(:) = Xreal(:, 2) + v;                                                  % Measurements and noise addition

%% Inputs
% u_1: Insulin [mU/min]
% u_2: Intravenous feed [mU/min]
% u_3: Enteral feed [mmol/min]

u = [u_1 u_2 u_3];                                                          % Inputs from simulink simulation
%% Particle filter parameters
Np = 100;                                                                 % Number of particles

i_filter = 1:5:size(Xreal, 1);                                              % Iterations for filter simulation
filter_size = size(i_filter, 2);
X_filter = zeros(7, filter_size);                                           % Memory allocation for filtered states

[Error, e_k_p] = noise_ARIMA(0, e_k_1);                                     % First noise generation with ARIMA assumption

X_filter(:, 1) = X_0;                                                 % First filtered state
X_filter(2, 1) = X_filter(2, 1) + Error;                                    % First measurement with system noise

% Tipos de errores a analizar
%RMSE = zeros(6, 1);

delta_X = zeros(6, filter_size);                                            % Error vector memory allocation
delta_X(:, 1) = Xreal(1, :) - X_filter([1:6], 1)';                          % First Error

step = i_filter(2) - i_filter(1);                                           % Integration time step
%% Particle filter implementation

tic;                                                                        % Time initialization
wk = repmat(1/Np, 1, Np);
for i = 2:size(i_filter, 2)   
    %wk = repmat(1/Np, [1, Np]);                                            % Particles weight initialization                                        % Weight initialization
    X_f_k_1 = X_filter(:, i-1);
    [X_filter(:, i), wk, e_k_p] = filtro(X_f_k_1, meas(i_filter(i)), ...    % Filtered state
                                         u(i_filter(i), :), wk, step, Np, e_k_p);
    delta_X(:, i) = Xreal(i_filter(i), :)' - X_filter([1:6], i);            % Error computation
end

RMSE_aux = sqrt(sum(delta_X .^ 2, 2) / filter_size);                        % E

tiempo = toc
RMSE = sum(RMSE_aux)
plot(tiempo_sim, Xreal(:, 1), i_filter, X_filter(1, :));