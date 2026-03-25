% CONTROLLED full order ABM
% Solve the optimal control problem in [t,t+nt*ht]
% Use controls ui(t+ht) to compute xi(t+ht)

clear all
close all
clc

% Initial conditions
load('data_d2_N10.mat');

% Parameters for the optimal control problem
gamma = 0.5;
itmax = 50; % maximum number of iteration
tol = 1e-6;

% Consensus parameter
tol_consensus = 1e-15;

% Time integration
t0 = 0;
tf = 1;
ht = 1e-2;
nt = length(t0:ht:tf);

% Interaction kernel
alpha = 1.6;
sig = @(y) 1./(1 + exp(-y));
phi = @(s) (1 - sig(alpha*(s - 1))) / (1 - sig(-alpha));
phi_prime = @(s) -(alpha/(1 - sig(-alpha))) .* sig(alpha*(s - 1)) .* (1 - sig(alpha*(s - 1)));

% Forward/Backward Sweep functions
Nag_cl = ones(N,1);
Hpfun = @(X,P) derivativeHp_op(X,P,N,Nag_cl,phi,gamma);
Hxfun = @(X,P) derivativeHx_op(X,P,N,Nag_cl,phi_prime,phi);

% Initialization
Xvec = [];
for i = 1:N
    Xvec = [Xvec, x{i}(:,1)];
end
X0 = Xvec';

mediax0 = zeros(N,1);
mediau0 = zeros(N,1);
for i = 1:N
    mediax0(i) = mean(x{i}(:,1));
end
xbar0 = zeros(d,1);
for i = 1 : N
    xbar0 = xbar0 + x{i}(:,1);
end
xbar0 = xbar0 / N;
Xt0 = 0;
for i = 1:N
    Xt0 = Xt0 + norm(x{i}(:,1) - xbar0)^2;
end
Xt0 = Xt0 / (2 * N^2);
mediax_all = mediax0;
mediau_all = mediau0;
Xt_all = Xt0;
t_all = 0;

% Iterative cycle
k = 0;
tic;
while Xt0 > tol_consensus
    k = k + 1;
    [Xout,Pout] = func_opt_control(N,d,nt,ht,X0,tol,itmax,Hpfun,Hxfun);
    X0 = [];
    xbar_val = 0;
    for i = 1 : N
        x{i}(:,k+1) = squeeze(Xout(2,i,:))';
        xbar_val = xbar_val + x{i}(:,k+1);
        Pend = squeeze(Pout(2,i,:))';
        u{i}(:,k) = -(N/(2*gamma))*Pend;
        mediax(i,k+1) = mean(x{i}(:,k+1));
        mediau(i,k) = mean(u{i}(:,k));
        X0 = [X0, x{i}(:,k+1)];
    end
    X0 = X0';
    xbar(:,k+1) = xbar_val./N;
    somma = 0;
    for i = 1 : N
        somma = somma + norm(x{i}(:,k+1) - xbar(:,k+1))^2;
    end
    Xt0 = somma/(2*N^2);
    X_t(k+1) = Xt0;
end
tf = t0 + ht*k;
tempo_totale = toc;
fprintf("Time: %.2f secondi\n", tempo_totale);

time = [0:ht:tf];

%% Plots
figure
semilogy(time,X_t,'LineWidth',2)
xlabel('t')
ylabel('X(t)')
title(['Consensus – \alpha = ', num2str(alpha)])
set(gca, 'XScale', 'log', 'FontSize', 12, 'FontWeight', 'B')
box on

figure
hold on
for i = 1 : N
    plot(time,mediax(i,:),'LineWidth',1)
    hold on
end
hold off
xlabel('t')
ylabel('Opinions')
title(['Controlled Opinion Dynamics – \alpha = ', num2str(alpha)])
set(gca, 'XScale', 'log', 'FontSize', 12, 'FontWeight', 'B')
box on


