function Hx = derivativeHx_op(X,P,N_tot,N_vec,phi_prime,phi)

K  = size(X,1);
r = size(X,2);
Xs = X;
Px = P;

Hx = zeros(K,r);
Xbar = zeros(1,r);
for i = 1 : K
    Xbar = Xbar + Xs(i,:);
end
Xbar = Xbar/K;

% average of the x-components
for i = 1 : K
    sumTermp = 0;
    for j = 1 : K
        if j == i, continue; end
        % Distance between "agent i" and "agent j"
        R_ij = Xs(i,:)-Xs(j,:);
        dist_ij = norm(R_ij);
        dotTerm = ((Px(j,:)-Px(i,:))*transpose((Xs(j,:)-Xs(i,:))))*(Xs(j,:)-Xs(i,:));
        val_ij = 0;
        if dist_ij > 1e-14
            val_ij = (phi_prime(dist_ij)/dist_ij)*dotTerm;
        end
        val_ij = val_ij + phi(dist_ij)*(Px(j,:)-Px(i,:));
        sumTermp = sumTermp + N_vec(j)*val_ij;
    end

    % Multiply by 1/N_tot and sum the " 2/K(x_i-xbar ) "
    Hx(i,:) = 1/N_tot*sumTermp + (2/K)*(Xs(i,:)-Xbar);
end
end
