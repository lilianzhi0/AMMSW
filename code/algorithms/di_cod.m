function [l_avg, l_ans, err_max, cost] = di_cod(X, Y, L_number, N, gap, flag_time, RX, RY)
    tic;
    [mx, n] = size(X);
    my = size(Y, 1);
    l_avg = 0;
    err_avg = 0;
    err_max = 0;
    l_sum = 0;
    l_ans = -1;
    cnt = 0;
    L = {};
    L_activate = {};

    % Initialise each layer
    for i = 1:L_number 
        L_activate{end + 1} = Block(mx, my, 2^(i+1), 0);
        L{end + 1} = {};
    end
    X_size = 0;
    Y_size = 0;
    len = 0;
    for t = 1:n
        if t >= N && mod(t, gap) == 0 && flag_time == 0
            [ll, epsilon] = query(L, t, N, X, Y);
            l_avg = l_avg + ll;
            err_avg = err_avg + epsilon;
            err_max = max(err_max,epsilon);
            cnt = cnt + 1;
            l_ans = max(l_ans, l_sum);
        end
        x = X(:, t);
        y = Y(:, t);
        for i = 1:L_number
            l = 2^i;
            if ~isempty(L{i}) && L{i}{1}.t_start <= t - N
                l_sum = l_sum - l;
                L{i}(1) =  [];
            end
            L_activate{i}.idx = L_activate{i}.idx + 1;
            L_activate{i}.X(:, (L_activate{i}.idx)) = x;
            L_activate{i}.Y(:, (L_activate{i}.idx)) = y;
            if L_activate{i}.idx >= 2 * l || t == n
                [L_activate{i}.X, L_activate{i}.Y, L_activate{i}.idx] = cod(1, L_activate{i}.X, L_activate{i}.Y, l);
            end 
            L_activate{i}.t_end = t;
        end
        X_size = X_size + norm(x, 2) ^ 2;
        Y_size = Y_size + norm(y, 2) ^ 2;
        
        if X_size > N * RX / 2^L_number || Y_size > N * RY / 2^L_number || t == n
            X_size = 0;
            Y_size = 0;
            len = len + 1;
            v = min(tail_zeros(len) + 1, L_number);
            for i = 1:v
                l = 2^i;
                [L_activate{i}.X, L_activate{i}.Y, L_activate{i}.idx] = cod(1, L_activate{i}.X, L_activate{i}.Y, l);
                block = Block(mx, my, l, t); 
                block.X = L_activate{i}.X(:, 1:l);
                block.Y = L_activate{i}.Y(:, 1:l);
                block.t_start = L_activate{i}.t_start;
                l_sum = l_sum + l;
                L{i}{end + 1} = block;
                L_activate{i}.clear(t); 
                L_activate{i}.t_start = t;
                L_activate{i}.t_end = t;
                L_activate{i}.X = zeros(L_activate{i}.mx, L_activate{i}.l);
                L_activate{i}.Y = zeros(L_activate{i}.my, L_activate{i}.l);
                L_activate{i}.idx = 0;
            end
        end
    end

    cost = toc;
    l_avg = l_avg / cnt;
    err_avg = err_avg / cnt;
end

function [ll, err] = query(L, t, N, X, Y)
    A = [];
    B = [];
    a = t;
    b = t - N;

    v = length(L);
    
    while v > 0 && isempty(L{v})
        v = v - 1;
    end
    
    st = 100000;
    en = -1;
    for i = v:-1:1
        if L{i}{1}.t_end <= a
            A = horzcat(L{i}{1}.X, A);
            B = horzcat(L{i}{1}.Y, B);
            a = L{i}{1}.t_start;
            st = min(st, a);
        end

        if length(L{i}) == 1
            b = L{i}{1}.t_end;
            en = max(en, b);
        elseif L{i}{end}.t_start >= b
            A = horzcat(A, L{i}{end}.X);
            B = horzcat(B, L{i}{end}.Y);
            b = L{i}{end}.t_end;
            en = max(en, b);
        end
    end

    XW = X(:, max(1, t - N + 1):t);
    YW = Y(:, max(1, t - N + 1):t);
    
    XYT = norm(XW, 'fro') * norm(YW, 'fro');
    err = norm_XYT_ABT(XW', YW', A', B') / XYT;
    [~, ll] = size(A);
end

function idx = tail_zeros(n)
    idx = 0;
    while bitget(n, 1) == 0
        idx = idx + 1;
        n = bitshift(n, -1);
    end
end
