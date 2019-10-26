% Team 20 - Avalanche Detection
% Nov 12th, algorithim demo
% Louis Rosenblum, Cayden Seiler, Khristian Jones

%% Sensor placement

s0 = [0 0];
s1 = [100 0];
s2 = [0 100];
s3 = [100 100];

%% Grid design

% data structure of all x,y locations for grid points
grid = cell(100,100);

for i = 1:100
    for j = 1:100
    grid{i,j} = [ (10*i-5) (10*j+995)];   
        
        
    end
end

%% Distance usage example
 
 dist1 = distance(s0,s1);
 
 dist1 = distance(s0, grid{30,80});



%% Avalanche condition generation

randx = randi(100,1,1);
randy = randi(100,1,1);

origin_point = {randx,randy}
origin = grid{randx, randy}

% Temp in celsius, -40 C to 10 C
tempc = randi([-40 10],1,1)
% tempc = tempk-273

% Universal gas constant
% r = 8.314;

% Adiabatic constant
% y = 1.4;

% Molecular mass for dry air
% m = .02895;

% Speed of sound in m/s
speed_of_sound = 331.3 * sqrt(1 + (tempc / 273.15))
% speed_of_sound = sqrt(y*r*tempk/m)




%% Calculate distance to sensors

d0 = distance(s0,origin) 
d1 = distance(s1,origin)
d2 = distance(s2,origin)
d3 = distance(s3,origin)

dist1 = distance(s0,grid{90,30});
delta1 = d1 - d0;
delta2 = d2 - d0;
delta3 = d3 - d0;

%% Signal Generation

figure();
t = 0:1/1024:0.25;


signal0 = cos(10*2*pi.*t);

plot(t,signal0), hold on

wavelength = speed_of_sound/10;
shift1 = delta1/wavelength
shift2 = delta2/wavelength;
shift3 = delta3/wavelength;
signal1 = cos(10*2*pi.*(t-shift1));
signal2 = cos(10*2*pi.*(t-shift2));
signal3 = cos(10*2*pi.*(t-shift3));

plot(t,signal1);
plot(t,signal2);
plot(t,signal3);

%% Algorithim definition




%% Plot 

figure();
% Sensors
scatter([0 0 100 100],[0 100 0 100],'filled'),xlim([-100 1100]),ylim([-100 2100]),hold on

% True origin
scatter([origin(1)],[origin(2)],'filled')

% Grid border
%plot([0 0 1000 1000 0],[1000 2000 2000 1000 1000],'g','Linewidth',2)

% Grid points
x1 = [];
y1 = [];

% One square filled to 100x100 resolution
for x = 1:10
    for y = 1:10
        z = grid{x,y};
        k1 = [(z(1) - 5) (z(1) +5) (z(1) +5) (z(1) -5) (z(1) -5)];
        k2 = [(z(2) + 5) (z(2) +5) (z(2) -5) (z(2) -5) (z(2) +5)];
        x1 = [x1 k1];
        y1 = [y1 k2];
    end
    plot(x1,y1,'b'), hold on;
    x1 = [];
    y1 = [];
end

% 10x10 resolution
for x = 1:10
    for y = 1:10
        z = grid{x*10,y*10};
        k1 = [(z(1) - 50) (z(1) +50) (z(1) +50) (z(1) -50) (z(1) -50)] - 45;
        k2 = [(z(2) + 50) (z(2) +50) (z(2) -50) (z(2) -50) (z(2) +50)] - 45;
        x1 = [x1 k1];
        y1 = [y1 k2];
    end
    plot(x1,y1,'b'),xlabel("m"),ylabel("m")
    x1 = [];
    y1 = [];
end



%% Gaussian noise


%% Distance function definition

function dist = distance(p1,p2)
    a = p2(1);
    b = p2(2);
    dist = sqrt(abs((p2(1) - p1(1))^2 + (p2(2)-p1(2))^2));
end