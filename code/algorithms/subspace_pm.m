function Z = subspace_pm(X, Y, m, q)
    dy = size(Y, 2);
    G = randn(dy, m);
    Z = X'*(Y*G);
    for i=1:q
        Z=X'*(Y*(Y'*(X*Z)));
        [Z,~]=qr(Z,0);
    end
end