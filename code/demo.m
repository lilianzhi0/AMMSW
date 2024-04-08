% This is a demo for experiments of DI-COD, DI-SCOD, EH-COD and EH-SCOD on dataset "APR".
% You can directly run "demo.m" to see the results.
% Due to the limit of file size, we only include dataset "APR" in supplementary materials.

setup;

filename = 'apr'; % set dataset
mm = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 30, 40, 50, 80, 100, 150, 200, 250, 300, 350, 400]; % set sketch size
ll = [6, 7, 8, 9, 10, 11]; % set level size
q = 5; % set iteration number of subspace power method

run_all_algorithms;
draw_all_figures;

% run_all_algorithms_time;
% draw_all_figures_time;
