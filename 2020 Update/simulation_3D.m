% 3D Simulation

%% Initialize program

clear all

%% Load 3D Elevation Map

latlim = [45.25532873 45.30078327]
longlim = [-111.4957325 -111.4048235]

[elevation, refvec] = dted("lone_peak.dt2",1,latlim,longlim)


%% Populate grid with X,Y, and Z coordinates

grid = cell(329,329);

for i = 1:329
    for k = 1:329
        x = (i-1)/328;
        y = (k-1)/328+2000;
        z = elevation(round(k*0.5),i);
        
        grid{i,k} = [x y z];
    end
    
end

%% Sensor placement

% Hardcoding sensor locations at (x,y) coordinates
s0 = [0 0 2885];
s1 = [100 0 2772];
s2 = [0 100 2648];
s3 = [100 100 2560];

%% Create countour plot

x = 0:1000/328:1000;
y = 2000:1000/164:3000;

contour(x,y,elevation,16)
colormap default

hold on

gscatter(0,0,'Sensor 0', 'b');
gscatter(0,100,'Sensor 1', 'r');
gscatter(100,0,'Sensor 2', 'y');
gscatter(100,100,'Sensor 3', 'm');

xlim([-100 1100]),ylim([-100 3100]);

hold off

    
%% Distance function definition

function dist = distance(p1,p2)
    a = p1(1);
    b = p1(2);
    c = p1(3);
    
    d = p2(1);
    e = p2(2);
    f = p3(3);

    dist = sqrt((d-a)^2+(e-b)^2+(f-c)^2);
end