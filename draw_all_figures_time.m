filepath = fullfile('output', filename);
L1 = readmatrix(fullfile(filepath, 'Sampling_time'));
L2 = readmatrix(fullfile(filepath, 'DI-COD_time'));
L3 = readmatrix(fullfile(filepath, 'DI-SCOD_time'));
L4 = readmatrix(fullfile(filepath, 'EH-COD_time'));
L5 = readmatrix(fullfile(filepath, 'EH-SCOD_time'));
L = {L1, L2, L3, L4, L5};

col = {'rx-', 'bo-', 'm^-', 'cs-', 'k*-', 'gv-', 'r<-'};

figure('Name', strcat(filepath, ': Average time vs. Average l'), 'NumberTitle', 'off');
for i = 1:(length(L))
    tmp = L{i}; avgl = tmp(1:end, 1); avgtime = tmp(1:end, 2); 
    plot(avgl, avgtime, col{i}, 'LineWidth', 1.2, 'MarkerSize', 8);
    hold on;
end
xlabel('Average l');
ylabel('Average time');
title(strcat(filename, ': Average time vs. Average l'));
lgd = legend('swor', 'DI', 'DI-sparse', 'EH', 'EH-sparse');

lgd.FontSize = 16;
lgd.Location = 'best';
set(gca, 'FontSize', 13);