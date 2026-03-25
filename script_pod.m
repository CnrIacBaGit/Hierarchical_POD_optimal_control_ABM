% Optimal control through model order reduction
% Reducing the dimension of each agent
% The subspace Psi_r onto which project the full model is updated every
% 100*ht steps (that corresponds to the time intervals [n,n+1]
clear all 
close all 
clc

% Initial conditions
load('data_d50_N50.mat')

% POD parameters
POD_interval = 100; % update frequency
tol_singval = 1e-3; % tolerance for the singular values
% Optimal control problem parameters
gamma = 0.5;
itmax = 50;
tol = 1e-6;
% Consensus parameter
tol_consensus = 1e-15;
% Time settings
t0 = 0; 
tf = 1; 
ht = 1e-2;
nt = length(t0:ht:tf);
% Communication kernel
alpha = 1.6;
sig = @(y) 1 ./ (1 + exp(-y));
phi = @(s) (1 - sig(alpha * (s - 1))) / (1 - sig(-alpha));
phi_prime = @(s) -(alpha / (1 - sig(-alpha))) .* sig(alpha * (s - 1)) .* (1 - sig(alpha * (s - 1)));

% Functions for the Forward/Backward Sweep
Nag_cl = ones(N,1);
Hpfun = @(X,P) derivativeHp_op(X,P,N,Nag_cl,phi,gamma);
Hxfun = @(X,P) derivativeHx_op(X,P,N,Nag_cl,phi_prime,phi);

xbar_val = 0;
for i = 1 : N
    mediax(i,1) = mean(x{i}(:,1));
    xbar_val = xbar_val + x{i}(:,1);
end
xbar(:,1) = xbar_val./N;

% Initialization
count = 0;
%% Solve the full model without control in [t0,tf]
tic;
count = count + 1;
fom_nocontrol;
for k = 1 : nt
    somma = 0;
    for i = 1 : N
        somma = somma + norm(x{i}(:,k) - xbar(:,k))^2;
    end
    Xt = somma/(2*N^2);
    X_t(k) = Xt;
end
%% First POD and projection
S = cat(2, x{:});
[Psi, Sigma, ~] = svd(S, 'econ');
r = sum(diag(Sigma) > tol_singval);
Psi_r = Psi(:, 1:r);
dim_r(count) = r;
for i = 1 : N
    xr{i} = Psi_r'*x{i}(:,end);
    xrfirst(i,1) = xr{i}(1);
end
%% Optimal control problem (in the reduced space)
% Initialization for the optimal control problem
Xr = [];
for i = 1 : N
    Xr = [Xr, xr{i}];
end
X0 = Xr';
[Xout,Pout] = func_opt_control(N,r,nt,ht,X0,tol,itmax,Hpfun,Hxfun);
for i = 1 : N
    Pend = squeeze(Pout(2,i,:))';
    uur = -(N/(2*gamma))*Pend;
    u{i} = Psi_r*uur';
    mediau(i,1) = mean(u{i});
    x{i} = x{i}(:,2:end);  % trim old state
end

%% Main loop
while Xt > tol_consensus
    Xt
    count = count + 1;
    % Step 1: advance system with current control
    fom_1step;
    somma = 0;
    for i = 1 : N
        somma = somma + norm(x{i}(:,end) - xbar(:,end))^2;
    end
    Xt = somma/(2*N^2);
    X_t = [X_t, Xt];
    % Step 2: update POD 
    if mod(count, POD_interval) == 0
        S = cat(2, x{:});
        [Psi, Sigma, ~] = svd(S, 'econ');
        r = sum(diag(Sigma) > tol_singval);
        Psi_r = Psi(:, 1:r);
    end
    dim_r(count) = r;
    for i = 1 : N
        xr{i} = Psi_r' * x{i}(:,end);
        xrfirst(i,count) = xr{i}(1);
    end
    Xr = [];
    for i = 1 : N
        Xr = [Xr, xr{i}];
    end
    X0 = Xr';
    % Step 3: optimal control in reduced space
    [Xout,Pout] = func_opt_control(N,r,nt,ht,X0,tol,itmax,Hpfun,Hxfun);
    for i = 1 : N
        Pend = squeeze(Pout(2,i,:))';
        uur = -(N/(2*gamma))*Pend;
        u{i} = Psi_r*uur';
        mediau(i,count) = mean(u{i});
        x{i} = x{i}(:,2:end);  % trim old state
    end
    T = tf + count*ht;
end
tempo_totale = toc;
fprintf("Time: %.2f secondi\n", tempo_totale);

time = [0:ht:T];
%% Plots
figure
semilogy(time(2:end),X_t,'LineWidth',2)
xlabel('t')
ylabel('X(t)')
title('Consensus parameter')
set(gca,'Fontsize',12,'FontWeight','B')

figure
hold on
for i = 1 : N
    plot(time(2:end),mediax(i,:),'LineWidth',1)
    hold on
end
hold off
xlabel('t')
title('Mean value for the position')
set(gca,'FontSize',12,'FontWeight','B')

figure
hold on
for i = 1 : N
    plot(time(nt+1:end),mediau(i,:),'LineWidth',1)
    hold on
end
hold off
xlabel('t')
title('Mean control effort')
set(gca,'FontSize',12,'FontWeight','B')

figure
plot(time(nt+1:end), dim_r, 'o');
xlabel('t') 
ylabel('r (reduced dimension)')
set(gca,'FontSize',12,'FontWeight','B')

figure
hold on
for i = 1 : N
    plot(time(2:end),xrfirst(i,:),'LineWidth',1)
    hold on
end
hold off
xlabel('t')
title('First coordinate reduced position')
set(gca,'FontSize',12,'FontWeight','B')