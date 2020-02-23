% Will put everything printed into the command window into a .txt file 
% to be used in calculating an average error.
diary('Diary.txt')
for i = 1:2
    simulation_3D()
end
diary('off')