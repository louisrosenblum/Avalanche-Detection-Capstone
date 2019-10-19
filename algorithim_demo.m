% Team 20 - Avalanche Detection
% Nov 12th, algorithim demo
% Louis Rosenblum, Cayden Seiler, Khristian Jones

%% Sensor initialization

% s0, reference sensor
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

%% Avalanche condition generation

randx = randi([1 100],1,1);
randy = randi([1 100],1,1);

origin = grid{randx, randy}

% Temp in kelvin, -40 F to 50 F
temp = randi([233 283],1,1)

% Universal gas constant
r = 8.314;

% Adiabatic constant
y = 1.4;

% Molecular mass for dry air
m = .02895;

% Speed of sound in m/s
speed_of_sound = sqrt(y*r*temp/m)

%% Calculate distance to sensors

d0 = distance(s0,origin)
d1 = distance(s1,origin)
d2 = distance(s2,origin)
d3 = distance(s3,origin)

delta1 = d1 - d0
delta2 = d2 - d0
delta3 = d3 - d0


%% Distance function definition

function dist = distance(p1,p2)
    dist = sqrt(abs((p2(1) - p1(1))^2 + (p2(2)-p1(2))^2));
end


