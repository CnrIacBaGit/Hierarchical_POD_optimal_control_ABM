% Solve the full model with explicit Euler scheme
% One step for the framework
xbar_val = 0;
for i = 1 : N
    somma = 0;
    for j = 1 : N
        dij = norm(x{i}(:,nt-1)-x{j}(:,nt-1));
        ck = phi(dij);
        somma = somma + ck.*(x{j}(:,nt-1)-x{i}(:,nt-1));
    end
    x{i}(:,nt) = x{i}(:,nt-1) + (ht/N).*somma + ht*u{i};
    mediax(i,nt+count-1) = mean(x{i}(:,nt));
    xbar_val = xbar_val + x{i}(:,nt);
end
xbar(:,nt) = xbar_val./N;

