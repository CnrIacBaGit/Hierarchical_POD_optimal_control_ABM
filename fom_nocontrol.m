% Solve the full model with explicit Euler scheme
for k = 1 : nt-1
    xbar_val = 0;
    for i = 1 : N
        somma = 0;
        for j = 1 : N
            dij = norm(x{i}(:,k)-x{j}(:,k));
            ck = phi(dij);
            somma = somma + ck.*(x{j}(:,k)-x{i}(:,k));
        end
        x{i}(:,k+1) = x{i}(:,k) + (ht/N).*somma;
        mediax(i,k+1) = mean(x{i}(:,k+1));
        xbar_val = xbar_val + x{i}(:,k+1);
    end
    xbar(:,k+1) = xbar_val./N;
end
