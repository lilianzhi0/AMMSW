function [l_avg, l_ans, err_max, cost] = eh_cod(X, Y, l, N, gap, flag_time)
    tic;
    
    [mx, n] = size(X);
    my = size(Y, 1);
    l_ans = -1;
    l_avg = l;
    err_avg = 0;
    err_max =0;
    b = l/2;
    q = 0;
    err = 0;
    maxsize = 0;
    L = {};
    lv1 = Level();
    L{end+1} = lv1;
    active = Block(mx, my, l, 0);
    
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
            l_ans = max(l_ans, sizeW);
            err_avg = err_avg + eW;
            err_max = max(err_max, eW);
        end
    
        x = X(:, t);
        y = Y(:, t);
        active = active.insert(x, y, t);
    
        if active.size >= l || t == n
            block = Block(mx, my, l, t);
            block.X = active.X(:, 1:active.idx);
            block.Y = active.Y(:, 1:active.idx);
            block.t_start = active.t_start;
            block.sizex = active.sizex;
            block.sizey = active.sizey;
            block.size = active.size;
            L{1}.blocks{end+1} = block;
            active = active.clear(t);
            for i = 1:length(L)
                if length(L{i}.blocks) > b
                    if i + 1 == length(L) + 1
                        new_level = Level();
                        L{end+1} = new_level;
                    end
                    if L{i}.blocks{1}.size > (2^(i))*l
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

