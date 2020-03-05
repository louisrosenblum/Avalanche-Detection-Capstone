% Will put everything printed into the command window into a .txt file 
% to be used in calculating an average error.
diary('log.txt')

result = []
tic

for go = 1:100
    go
    simulation_3D();
    result(go) = max_z;
end
diary('off')

toc

