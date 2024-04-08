% the runtime of algorithms
load(filename);
X = X';
Y = Y';
norm_of_XYT = norm_XYT(X', Y');

[mx, n] = size(X);
my = size(Y, 2);
N = 10000; % the size of window 
gap = 10000;
flag_time = 1;

l_avg = zeros(5, max(length(ll), length(mm)));
cost = zeros(5, max(length(ll), length(mm)));

for i=1: length(ll)
    L = ll(i);
    fprintf('data=%s N=%d level = %d\n', filename, N, L);
    l_avg(1, i) = L;
    l_avg(2, i) = L;
    
    [~, ~, ~, cost(1, i)] = di_cod(X, Y, L, N, gap, flag_time, RX, RY);
    fprintf('DI-COD: L = %d, time = %f\n', L, cost(1, i));

    [~, ~, ~, cost(2, i)] = di_scod(X, Y, L, N, gap, q, flag_time, RX, RY);
    fprintf('DI-SCOD: L = %d, time = %f\n', L, cost(2, i));
end

for i=1: length(mm)
    m = mm(i);
    fprintf('data=%s N=%d m = %d\n', filename, N, m);
    l_avg(3, i) = m;
    l_avg(4, i) = m;
    l_avg(5, i) = m;

    [~, ~, ~, cost(3, i)] = eh_cod(X, Y, m, N, gap, flag_time);
    fprintf('EH-COD: l = %d, time = %f\n', m, cost(3, i));

    [~, ~, ~, cost(4, i)] = eh_scod(X, Y, m, N, gap, q, flag_time);
    fprintf('EH-SCOD: l = %d, time = %f\n', m, cost(4, i));

    [~, ~, ~, cost(5, i)] = sampling(X, Y, m, N, gap, flag_time);
    fprintf('Sampling: l = %d, time = %f\n', m, cost(5, i));
    
end

fn = {'DI-COD', 'DI-SCOD', 'EH-COD', 'EH-SCOD', 'Sampling'};

for i = 1:5
    tmp = 0;
    if i <= 2
        tmp = length(ll);
    else
        tmp = length(mm);
    end
    load = fullfile('output', filename);
    fileID = fopen(fullfile(load, sprintf('%s_time.txt', fn{i})), 'w');
    fprintf(fileID, 'l, time\n');
    for j = 1:tmp
        fprintf(fileID, '%d, %f\n', floor(l_avg(i, j)), cost(i, j));
    end
    fclose(fileID);
end

