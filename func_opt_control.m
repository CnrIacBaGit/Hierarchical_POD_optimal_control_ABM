function [Xout,Pout] = func_opt_control(N,d,nt,ht,X0,tol,itmax,Hpfun,Hxfun)

Xout = zeros(nt, N, d);
Pout = zeros(nt, N, d);
PT = zeros(N,d);
par = 1;
for m = 1 : nt
    Xout(m,:,:) = reshape(X0,1,N,d);
    Pout(m,:,:) = reshape(PT,1,N,d);
end
Xold = Xout;
Pold = Pout;
stop = false;
it = 0;
while ~stop
    it = it + 1;
    XX = X0;
    Xout(1,:,:) = reshape(XX,1,N,d);
    for n = 2 : nt
        P_current = squeeze(Pout(n,:,:));
        if d == 1
            P_current = P_current';
        end
        Hp_val = Hpfun(XX, P_current);
        XX = XX + ht * Hp_val;
        Xout(n,:,:) = XX;
    end
    P = PT;
    Pout(nt,:,:) = reshape(P,1,N,d);
    for n = nt-1:-1:1
        X_current = squeeze(Xout(n,:,:));
        if d == 1
            X_current = X_current';
        end
        Hx_val = Hxfun(X_current, P);
        P = par * (P + ht * Hx_val);
        Pout(n,:,:) = P;
    end
    rel_err = max(norm(Xout(:) - Xold(:)) / max(norm(Xout(:)), 1e-14), norm(Pout(:) - Pold(:)) / max(norm(Pout(:)), 1e-14));
    stop = (rel_err < tol) || (it >= itmax);
    Xold = Xout;
    Pold = Pout;
end