function sigma_1 = norm_XYT(X, Y, T)
    if nargin < 3
        T = 50;
    end
        
    v = randn(size(Y, 2), 1);
    v = v / norm(v);
    
    for i=1: T
        v = Y'*(X*(X'*(Y*v))); 
        v = v / norm(v);
    end
    sigma_1 = norm(X'*(Y*v));
end