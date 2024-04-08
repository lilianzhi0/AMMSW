filepath = fullfile('output', filename);
L1 = readmatrix(fullfile(filepath, 'Sampling'));
L2 = readmatrix(fullfile(filepath, 'DI-COD'));
L3 = readmatrix(fullfile(filepath, 'DI-SCOD'));
L4 = readmatrix(fullfile(filepath, 'EH-COD'));
L5 = readmatrix(fullfile(filepath, 'EH-SCOD'));
L = {L1, L2, L3, L4, L5};

% 1 sketch size 2 framework size

flag = 0; 
col = {'rx-', 'bo-', 'm^-', 'cs-', 'k*-', 'gv-', 'r<-'};
    
if flag == 1 || flag == 0
    figure('Name', strcat(filepath, ': Average Error vs. Average l'), 'NumberTitle', 'off');
    for i = 1:(length(L))
        tmp = L{i}; avgl = tmp(1:end, 1); avgError = tmp(1:end, 3); 
        plot(avgl, avgError, col{i}, 'LineWidth', 1.2, 'MarkerSize', 8);
        hold on;
    end
    xlabel('Average l');
    ylabel('Average Error');
    title(strcat(filename, ': Average Error vs. Average l'));
    lgd = legend('swor', 'DI', 'DI-sparse', 'EH', 'EH-sparse');
elseif flag == 2 || flag == 0
    figure('Name', strcat(filepath, ': Average Error vs. Max sketch size (columns)'), 'NumberTitle', 'off');
    for i = 1:(length(L))
        tmp = L{i}; avgl = tmp(1:end, 2); avgError = tmp(1:end, 3); 
        plot(avgl, avgError, col{i}, 'LineWidth', 1.2, 'MarkerSize', 8);
        hold on;
    end
    xlabel('Max sketch size (columns)');
    ylabel('Average Error');
    title(strcat(filename, ': Average Error vs. Max sketch size (columns)'));
    lgd = legend('swor', 'DI', 'DI-sparse', 'EH', 'EH-sparse');
end

lgd.FontSize = 16;
lgd.Location = 'best';
set(gca, 'FontSize', 13);