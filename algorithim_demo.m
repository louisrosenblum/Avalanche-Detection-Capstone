% Team 20 - Avalanche Detection
% Nov 12th, algorithim demo
% Louis Rosenblum, Cayden Seiler, Khristian Jones

%% Sensor initialization

s0 = [0 0];
s1 = [100 0];
s2 = [0 100];
s3 = [100 100];

%% Grid intialization

% Initialize grid
grid = cell(100,100);

% Populate grid points with x,y coordinates (center of each grid point)
for i = 1:100
    for j = 1:100
    grid{i,j} = [ (10*i-5) (10*j+995)];     
    end
end

%% Distance usage example

dist1 = distance(s0,grid{90,30})


%% Distance function definition

function dist = distance(p1,p2)
    dist = sqrt(abs((p2(1) - p1(1))^2 + (p2(2)-p1(2))^2));
end


