function [l_avg, l_ans, err_max, cost] = di_scod(X, Y, L_number, N, gap, q, flag_time, RX, RY)
    tic;
    [mx, n] = size(X);
    my = size(Y, 1);
    
    l_avg = 0;
    err_avg = 0;
    err_max = 0;
    l_ans = -1;
    l_sum = 0;
    cnt = 0;

    L = {};
    L_activate = {};

    [col_idx_X, row_idx_X, value_X] = find(X);
    [col_idx_Y, row_idx_Y, value_Y] = find(Y);

    nnz_row_X = full(sum(X~=0,1)); % nnz of each row of X 
    nnz_row_Y = full(sum(Y~=0,1)); % nnz of each row of Y 

    % Initialise each layer
    for i = 1:L_number
        L_activate{end + 1} = Buffer(mx, my);
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
            err_max = max(err_max, epsilon);
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

            L_activate{i}.X_nnz = L_activate{i}.X_nnz + nnz_row_X(t);
            L_activate{i}.Y_nnz = L_activate{i}.Y_nnz + nnz_row_Y(t);
            L_activate{i}.columns = L_activate{i}.columns + 1;
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

                % sparse matrix
                L_activate{i}.X_buffer = sparse(row_idx_X(L_activate{i}.X_entry_start: L_activate{i}.X_entry_start+L_activate{i}.X_nnz-1)-row_idx_X(L_activate{i}.X_entry_start)+1,...
                    col_idx_X(L_activate{i}.X_entry_start: L_activate{i}.X_entry_start+L_activate{i}.X_nnz-1),...
                    value_X(L_activate{i}.X_entry_start: L_activate{i}.X_entry_start+L_activate{i}.X_nnz-1),...
                    L_activate{i}.columns, L_activate{i}.mx);
                L_activate{i}.Y_buffer = sparse(row_idx_Y(L_activate{i}.Y_entry_start: L_activate{i}.Y_entry_start+L_activate{i}.Y_nnz-1)-row_idx_Y(L_activate{i}.Y_entry_start)+1,...
                    col_idx_Y(L_activate{i}.Y_entry_start: L_activate{i}.Y_entry_start+L_activate{i}.Y_nnz-1),...
                    value_Y(L_activate{i}.Y_entry_start: L_activate{i}.Y_entry_start+L_activate{i}.Y_nnz-1),...
                    L_activate{i}.columns, L_activate{i}.my);

                % Construct tX and tY by Subspace Power Itearation 
                Z = subspace_pm(L_activate{i}.X_buffer, L_activate{i}.Y_buffer, 2 * l, q);
                [tU, tS, tV]=svd((Z' * L_activate{i}.X_buffer') * L_activate{i}.Y_buffer, 'econ');
                ts=diag(tS);
                tX = diag(ts.^0.5) * tU' * Z';
                tY = diag(ts.^0.5) * tV'; 
                
                block = Block(mx, my, l, t); 
                [block.X, block.Y, ~] = cod(1, tX', tY', l);
                block.t_start = L_activate{i}.t_start;
                L{i}{end + 1} = block;
                l_sum = l_sum + l;
                L_activate{i}.reset(t);
                L_activate{i}.X_entry_start = L_activate{i}.X_entry_start + L_activate{i}.X_nnz;
                L_activate{i}.Y_entry_start = L_activate{i}.Y_entry_start + L_activate{i}.Y_nnz;
                L_activate{i}.X_nnz = 0;
                L_activate{i}.Y_nnz = 0;
                L_activate{i}.columns = 0;
                L_activate{i}.t_start = t;
                L_activate{i}.t_end = t;
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
