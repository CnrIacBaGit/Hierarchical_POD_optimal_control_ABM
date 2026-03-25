% Optimal control 
% Reduction in the number of agents by clustering

clear all
close all
clc

% Initial conditions
load('data_d50_N50.mat')

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
% Hpfun = @(X,P) derivativeHp_op(X,P,N,Nag_cl,phi,gamma);
% Hxfun = @(X,P) derivativeHx_op(X,P,N,Nag_cl,phi_prime,phi);

% Initialization
count = 0;
for i = 1 : N
    mediax(i,1) = mean(x{i}(:,1));
end
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
% Clustering by using DBSCAN
X = [];
for i = 1 : N
    X = [X, x{i}(:,end)];
end
epsilon = norm(X)/N;
index = dbscan(X',epsilon,1);
K = max(index);
cluster(count) = K; % save the number of clusters at each iteration
if K <= 1
    K = 2;
    [index,C] = kmeans(X',K);
end
% Find agents in each cluster and form the new variables xhat, vhat
xhat = [];
Nag_cl = [];
for i = 1 : K
    sommax = 0;
    indici = find(index==i);
    Ni = length(indici); % How many agents in cluster i?
    Nag_cl(i) = Ni;
    for j = 1 : Ni
        jj = indici(j);
        sommax = sommax + x{jj}(:,end);
    end
    xhat{i}(:,1) = sommax/Ni;
end
Xhat = [];
for i = 1 : K
    Xhat = [Xhat, xhat{i}];
end
X0 = Xhat';
Hpfun = @(X,P) derivativeHp_op(X,P,N,Nag_cl,phi,gamma);
Hxfun = @(X,P) derivativeHx_op(X,P,N,Nag_cl,phi_prime,phi);
[Xout,Pout] = func_opt_control(K,d,nt,ht,X0,tol,itmax,Hpfun,Hxfun);
% Define u{i} the control for each agent i
% u{i} = uhat{k} with i in Ik (cluster k)
for i = 1 : N
    idx = index(i);
    Pend = squeeze(Pout(2,idx,:))';
    uur = -(K/(2*gamma))*Pend;
    u{i} = uur';
    mediau(i,1) = mean(u{i});
    x{i} = x{i}(:,2:end);
end

while Xt > tol_consensus 
    Xt
    count = count + 1;
    % Step 1: advance system with current control
    fom_1step;
    somma = 0;
    for i = 1 : N
        somma = somma + norm(x{i}(:,end)-xbar(:,end))^2; % Compute Gamma for consensus
    end
    Xt = somma./(2*N^2);
    X_t = [X_t, Xt];
    % Step 2: Clustering by using DBSCAN
    X = [];
    for i = 1 : N
        X = [X, x{i}(:,end)];
    end
    index = dbscan(X',epsilon,1);
    K = max(index);
    cluster(count) = K; % save the number of clusters at each iteration
    if K <= 1
        K = 2;
        [index,C] = kmeans(X',K);
    end
    % Find agents in each cluster and form the new variables xhat, vhat
    xhat = [];
    Nag_cl = [];
    for i = 1 : K
        sommax = 0;
        indici = find(index==i);
        Ni = length(indici); % How many agents in cluster i?
        Nag_cl(i) = Ni;
        for j = 1 : Ni
            jj = indici(j);
            sommax = sommax + x{jj}(:,end);
        end
        xhat{i}(:,1) = sommax/Ni;
    end
    % Step 3: Solve the optimal control problem
    Xhat = [];
    for i = 1 : K
        Xhat = [Xhat, xhat{i}];
    end
    X0 = Xhat';
    Hpfun = @(X,P) derivativeHp_op(X,P,N,Nag_cl,phi,gamma);
    Hxfun = @(X,P) derivativeHx_op(X,P,N,Nag_cl,phi_prime,phi);
    [Xout,Pout] = func_opt_control(K,d,nt,ht,X0,tol,itmax,Hpfun,Hxfun);
    % Define u{i} the control in the full space for each agent i
    % u{i} = uhat{k} with i in Ik (cluster k)
    for i = 1 : N
        idx = index(i);
        Pend = squeeze(Pout(2,idx,:))';
        uur = -(K/(2*gamma))*Pend;
        u{i} = uur';
        mediau(i,count) = mean(u{i});
        x{i} = x{i}(:,2:end);
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
plot(time(nt+1:end),cluster,'o')
xlabel('t')
ylabel('Number of clusters')
set(gca,'Fontsize',12,'FontWeight','B')