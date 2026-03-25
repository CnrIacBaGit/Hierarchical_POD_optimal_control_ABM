% Full order ABM without control

clear all
close all
clc

% Load initial data
load('data_d2_N10.mat')

alpha = 0.1;     % HK parameter
ht = 1e-3;       % time step
tol = 1e-19;     % consensus threshold
t0 = 0;          % initial time
tmax = 30;     % final time
itmax = length(t0:ht:tmax);

% Interaction kernel
sig = @(y) 1./(1 + exp(-y));
phi = @(s) (1 - sig(alpha*(s - 1))) / (1 - sig(-alpha));

% Initialization
k = 1;
mediax = zeros(N,1);
X_t = [];

% Compute initial xbar
xbar_val = 0;
for i = 1 : N
    xbar_val = xbar_val + x{i}(:,1);
end
xbar(:,1) = xbar_val / N;

% Compute inital Xt
somma = 0;
for i = 1:N
    somma = somma + norm(x{i}(:,1) - xbar(:,1))^2;
    mediax(i,1) = mean(x{i}(:,1));
end
X_t(1) = somma / (2*N^2);

% Evolution
while X_t(end) > tol &&  k < itmax
    k = k + 1;
    xbar_val = 0;
    somma = 0;
    for i = 1 : N
        interazione = 0;
        for j = 1 : N
            if j ~= i
                dij = norm(x{j}(:,k-1) - x{i}(:,k-1));
                interazione = interazione + phi(dij) * (x{j}(:,k-1) - x{i}(:,k-1));
            end
        end
        x{i}(:,k) = x{i}(:,k-1) + (ht/N) * interazione;
        xbar_val = xbar_val + x{i}(:,k);
        mediax(i,k) = mean(x{i}(:,k));
    end
    xbar(:,k) = xbar_val / N;
    for i = 1 : N
        somma = somma + norm(x{i}(:,k) - xbar(:,k))^2;
    end
    X_t(k) = somma / (2*N^2);
end
tf = t0 + ht*k;
tspan = t0:ht:tf-ht;

%% Plot consensus
figure
semilogy(tspan, X_t, 'LineWidth', 2)
xlabel('t')
ylabel('X(t)')
title(['Consensus – \alpha = ', num2str(alpha)])
set(gca, 'XScale', 'log', 'FontSize', 12, 'FontWeight', 'B')
box on

%% Plot opinions
figure
hold on
for i = 1:N
    plot(tspan, mediax(i,:), 'LineWidth', 1)
end
hold off
xlabel('t')
ylabel('Opinions')
title(['Opinion Dynamics – \alpha = ', num2str(alpha)])
set(gca, 'XScale', 'log', 'FontSize', 16, 'FontWeight', 'B')
box on