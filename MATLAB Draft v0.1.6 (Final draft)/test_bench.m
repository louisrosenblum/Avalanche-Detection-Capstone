% Will put everything printed into the command window into a .txt file 
% to be used in calculating an average error.
diary('log.txt')

fail_ = []
pass_ = []

tic

for go = 1:10
    go
    simulation_3D();
    fail_(go) = fail;
    pass_(go) = pass;
end
diary('off')

toc

average_fail = mean(fail_)

average_pass = mean(pass_)