function [err,sizeW] = eh_query(L, l, t, X, Y, N, mx, my)
    qb = Block(mx, my, l, 0);
    qb.X = zeros(mx, 0);
    qb.Y = zeros(my, 0);
    sizeW = 0;
    for i = 1:length(L)
        for j = 1:length(L{i}.blocks)
            sizeW = sizeW + size(L{i}.blocks{j}.X,2);
            [A, B] = merge(qb, L{i}.blocks{j}, l);
            qb.X = A;
            qb.Y = B;
        end
    end
    XW = X(:, max(1, t-N+1):t);
    YW = Y(:, max(1, t-N+1):t);
    XFYF = norm(XW, 'fro') * norm(YW, 'fro');
    err = norm_XYT_ABT(XW', YW', qb.X', qb.Y') / XFYF; 
end

function [A, B] = merge(block1, block2, l)
    X = horzcat(block1.X, block2.X);
    Y = horzcat(block1.Y, block2.Y);
    if size(X, 2) <= l
        A = X;
        B = Y;
    else
        [A, B] = cod(2, X, Y, l); 
    end
end