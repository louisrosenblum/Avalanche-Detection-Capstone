% Will put everything printed into the command window into a .txt file 
% to be used in calculating an average error.
diary('log.txt')

results = []

tic

for go = 1:100
    go
    simulation_3D();
    results(go) = error;
end
diary('off')

toc

average_error = mean(results)

std(results)