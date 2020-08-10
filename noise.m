clear all; clc;

mu_ARIMA = 0;                                                               % Measurement noise mean with ARIMA assumption
sigma_ARIMA = sqrt(0.1);                                                    % Measurement noise variance with ARIMA assumption
e_k_1 = randn().*sigma_ARIMA + mu_ARIMA;                                    % First e value in ARIMA assumption

N = 10000;

for i = 1:N
    [error(i), e_k_1] = noise_ARIMA(i, e_k_1);
end

norm_error = error/sum(error);

pd = fitdist(norm_error','Normal')

hold on;
histogram(norm_error);
x = [-0.1:.01:0.1];
norm = normpdf(x, pd.mu, 0.25*pd.sigma);
plot(x,norm);
hold off;


