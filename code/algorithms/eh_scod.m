function [l_avg, l_ans, err_max, cost] = eh_scod(X, Y, l, N, gap, q, flag_time)
    tic;

    [col_idx_X, row_idx_X, value_X] = find(X);
    [col_idx_Y, row_idx_Y, value_Y] = find(Y);
    nnz_row_X = full(sum(X~=0,1)); % nnz of each row of X 
    nnz_row_Y = full(sum(Y~=0,1)); % nnz of each row of Y 
    
    [mx, n] = size(X);
    my = size(Y, 1);
    l_avg = l;
    l_ans = -1;
    err_avg = 0;
    err_max = 0;

    X_entry_start = 1;
    Y_entry_start = 1;
    buffer_X_nnz = 0;
    buffer_Y_nnz = 0;
    buffer_row_number = 0;
    buffer_start = 0;
    buffer_X_size = 0;
    buffer_Y_size = 0;
    buffer_size = 0;
    b = l/2;
    err = 0;
    maxsize = 0;
    L = {};
    lv1 = Level();
    L{end+1} = lv1;
    
    for t = 1:n
    
        % Pop expired blocks
        if  ~isempty(L{end}.blocks) && L{end}.blocks{1}.t_end < t - N
            L{end}.blocks(1) = [];
            if isempty(L{end}.blocks)
                L(end) = [];
            end
        end
    
        % Query
        if t >= N && mod(t, gap) == 0 && flag_time == 0
            q = q + 1;
            [eW, sizeW] = eh_query(L, l, t, X, Y, N, mx, my); 
            err_avg = err_avg + eW;
            err_max = max(err_max, eW);
            l_ans = max(l_ans, sizeW);
        end
    
        x = X(:, t);
        y = Y(:, t);
        
        buffer_X_nnz = buffer_X_nnz + nnz_row_X(t);
        buffer_Y_nnz = buffer_Y_nnz + nnz_row_Y(t);
        buffer_row_number = buffer_row_number + 1;
        buffer_X_size = buffer_X_size + norm(x,2)^2;
        buffer_Y_size = buffer_Y_size + norm(y,2)^2;
        buffer_size = sqrt(buffer_X_size) * sqrt(buffer_Y_size);
    
        if buffer_size >= 10*l || t == n
            X_buffer = sparse(row_idx_X(X_entry_start: X_entry_start+buffer_X_nnz-1)-row_idx_X(X_entry_start)+1,...
                    col_idx_X(X_entry_start: X_entry_start+buffer_X_nnz-1),...
                    value_X(X_entry_start: X_entry_start+buffer_X_nnz-1),...
                    buffer_row_number, mx);
            Y_buffer = sparse(row_idx_Y(Y_entry_start: Y_entry_start+buffer_Y_nnz-1)-row_idx_Y(Y_entry_start)+1,...
                    col_idx_Y(Y_entry_start: Y_entry_start+buffer_Y_nnz-1),...
                    value_Y(Y_entry_start: Y_entry_start+buffer_Y_nnz-1),...
                    buffer_row_number, my);
            % Construct tX and tY by Subspace Power Itearation
            if buffer_row_number > 4*l
                m = l;
            else
                m = ceil(buffer_row_number/4);
            end
            if m == buffer_row_number
                tX = full(X_buffer);
                tY = full(Y_buffer);
            else
                Z = subspace_pm(X_buffer, Y_buffer, m, q);
                [tU, tS, tV]=svd((Z' * X_buffer') * Y_buffer, 'econ');
                ts=diag(tS);
                tX = diag(ts.^0.5) * tU' * Z';
                tY = diag(ts.^0.5) * tV';
            end
            block = Block(mx, my, l, t);
            block.X = tX';
            block.Y = tY';
            block.t_start = buffer_start;
            block.sizex = buffer_X_size;
            block.sizey = buffer_Y_size;
            block.size = buffer_size;
            L{1}.blocks{end+1} = block;
            
            % reset buffer        
            X_entry_start = X_entry_start + buffer_X_nnz;
            Y_entry_start = Y_entry_start + buffer_Y_nnz;
            buffer_X_nnz = 0;
            buffer_Y_nnz = 0;
            buffer_row_number = 0;
            buffer_X_size = 0;
            buffer_Y_size = 0;
            buffer_size = 0;
            buffer_start = t;

            for i = 1:length(L)
                if length(L{i}.blocks) > b
                    if i + 1 == length(L) + 1
                        new_level = Level();
                        L{end+1} = new_level;
                    end
                    if L{i}.blocks{1}.size > (2^(i))*l*10
                        block0 = L{i}.blocks{1};
                        L{i}.blocks(1) = [];
                        L{i + 1}.blocks{end+1} = block0;
                    else
                        new_block = Block(mx, my, l, t);
                        new_block.t_start = L{i}.blocks{1}.t_start;
                        new_block.t_end = L{i}.blocks{2}.t_end;
                        new_block.sizex = L{i}.blocks{1}.sizex + L{i}.blocks{2}.sizex;
                        new_block.sizey = L{i}.blocks{1}.sizey + L{i}.blocks{2}.sizey;
                        new_block.size = sqrt(new_block.sizex) * sqrt(new_block.sizey);
                        XX = horzcat(L{i}.blocks{1}.X, L{i}.blocks{2}.X);
                        YY = horzcat(L{i}.blocks{1}.Y, L{i}.blocks{2}.Y);
                        if size(XX, 2) <= l
                            new_block.X = XX;
                            new_block.Y = YY;
                        else
                            [new_block.X, new_block.Y] = cod(2, XX, YY, l); 
                        end
                        L{i}.blocks(1:2) = [];
                        L{i + 1}.blocks{end+1} = new_block;
                        
                    end
                end
            end
        end
    end
    err_avg = err_avg / q;
    cost = toc;
end