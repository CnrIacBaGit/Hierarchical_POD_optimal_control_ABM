function Hp = derivativeHp_op(X,P,N_tot,N_vec,phi,gamma)

K  = size(X,1);
r = size(X,2);

Xs = X;
Px = P;

Hp = zeros(K,r);
u = -(K/(2*gamma)).*Px;
for i = 1 : K
    sumTerm = 0;
    for j = 1 : K
        if j == i, continue; end
        % Distance between "agent i" and "agent j"
        R_ij = Xs(i,:) - Xs(j,:);
        dist_ij = norm(R_ij);
        val_ij = phi(dist_ij)*(Xs(j,:) - Xs(i,:));
        sumTerm = sumTerm + N_vec(j)*val_ij;
    end
    Hp(i,:) = (1/N_tot)*sumTerm + u(i,:);
end
end
