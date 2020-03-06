% Simulink Testbench Script
% Flat Earth Inc
% Louis Rosenblum and Khristian Jones
% 03/06/2020

% Khristian see section titled 'Function call' for phase delays, signal
% data stored in signal0 - signal7
%% Load 3D Elevation Map

latlim = [45.25532873 45.30078327];
longlim = [-111.4957325 -111.4048235];

[elevation, refvec] = dted("lone_peak.dt2",1,latlim,longlim);

%% Populate grid with X,Y, and Z coordinates

grid = cell(150,150);

for i = 1:150
    for k = 1:150
        x = 1000*(i-1)/149;
        y = 1000*(k-1)/149;
        z = elevation(round(k*0.5*329/150),round(i*329/150));
        grid{i,k} = [x y z];

    end

end

%% Sensor placement

% Hardcoding sensor locations at (x,y) coordinates
s0 = [0 -1000 2596];
s1 = [100 -1000 2749];

s2 = [400 -1000 2671];
s3 = [500 -1000 2624];

s4 = [0 -1100 2613];
s5 = [100 -1100 2807];

s6 = [400 -1100 2726];
s7 = [500 -1100 2728];

%% Signal condition generation

% Two random intergers from 1-100 for origin out of possible grid indexes
randx = randi(150,1,1);
randy = randi(150,1,1);

origin = grid{randx, randy};

% Generate random temp in celsius, -40 C to 10 C
tempc = randi([-40 10],1,1)

% Calculate speed of sound in m/s
speed_of_sound = 331.3 * sqrt(1 + (tempc / 273.15))

%% Calculate distance to sensors

d0 = distance(origin,s0);
d1 = distance(origin,s1);
d2 = distance(origin,s2);
d3 = distance(origin,s3);

d4 = distance(origin,s4);
d5 = distance(origin,s5);
d6 = distance(origin,s6);
d7 = distance(origin,s7);

% Calculate difference in distance from sensors 1-3 to reference sensor 0
delta1 = d1 - d0;
delta2 = d2 - d0;
delta3 = d3 - d0;
delta4 = d4 - d0;
delta5 = d5 - d0;
delta6 = d6 - d0;
delta7 = d7 - d0;

% Calculate amplitude decay over each distance based on energy distributed
% over surface area of a sphere
decay0 = 100000000/(4*pi*d0^2);
decay1 = 100000000/(4*pi*d1^2);
decay2 = 100000000/(4*pi*d2^2);
decay3 = 100000000/(4*pi*d3^2);
decay4 = 100000000/(4*pi*d4^2);
decay5 = 100000000/(4*pi*d5^2);
decay6 = 100000000/(4*pi*d6^2);
decay7 = 100000000/(4*pi*d7^2);

%% Signal Generation

% Sampling rate
f = 1000

% Time vector
t = -4:1/f:4;

wave = rand(1,length(t));

% Generate signal hitting the reference sensor
signal0 = decay0 .* sinc(10*2*pi.*t) .* heaviside(t).*wave;

% Calculate wavelength
wavelength = speed_of_sound/10;

% Shift other signals to match distance travelled to each sensor
shift1 = delta1/wavelength;
shift2 = delta2/wavelength;
shift3 = delta3/wavelength;
shift4 = delta4/wavelength;
shift5 = delta5/wavelength;
shift6 = delta6/wavelength;
shift7 = delta7/wavelength;

% Generate signals received by each sensor
signal1 = decay1 .* sinc(10*2*pi.*(t-shift1/10)).* heaviside(t-shift1/10).*wave;
signal2 = decay2 .* sinc(10*2*pi.*(t-shift2/10)).* heaviside(t-shift2/10).*wave;
signal3 = decay3 .* sinc(10*2*pi.*(t-shift3/10)).* heaviside(t-shift3/10).*wave;

signal4 = decay4 .* sinc(10*2*pi.*(t-shift4/10)).* heaviside(t-shift4/10).*wave;
signal5 = decay5 .* sinc(10*2*pi.*(t-shift5/10)).* heaviside(t-shift5/10).*wave;
signal6 = decay6 .* sinc(10*2*pi.*(t-shift6/10)).* heaviside(t-shift6/10).*wave;
signal7 = decay7 .* sinc(10*2*pi.*(t-shift7/10)).* heaviside(t-shift7/10).*wave;

% Power factor
pf = 0;
snr = round(rand(1)*24 + 6)

% Add independent gaussian noise to each signal
signal0 = awgn(signal0,snr,pf);
signal1 = awgn(signal1,snr,pf);
signal2 = awgn(signal2,snr,pf);
signal3 = awgn(signal3,snr,pf);

signal4 = awgn(signal4,snr,pf);
signal5 = awgn(signal5,snr,pf);
signal6 = awgn(signal6,snr,pf);
signal7 = awgn(signal7,snr,pf);

%% Function call

% f - sample rate
% s0 - reference sensor
% s1 - alternate sensor
% s2 - grid structure
% x - horizontal index
% y - vertical index

for x = 1:150
    for y = 1:150
        phase_delay = delay(f,s0,s1,grid,x,y);
    end   
end

%% Fetch delays function

function shift = delay(f,s0,s1,grid,i,k)
            
            % Calculate distance from current grid point to each sensor
            distance0 = distance(s0,grid{i,k});
            distance1 = distance(s1,grid{i,k});
            
            % Determine difference in distance to reach sensor 1-3 compared
            % to reference sensor 0
            delta_1 = distance1 - distance0;

            
            % Calculate wavelength from speed of sound
            wave_length = speed/10;
            
            % Calculate phase shifts from wavelength
            shift = delta_1/wave_length*f;
            
end

%% Distance function definition

function dist = distance(p1,p2)
    a = p1(1);
    b = p1(2);
    c = p1(3);
    
    d = p2(1);
    e = p2(2);
    f = p2(3);

    dist = sqrt((d-a)^2+(e-b)^2+(f-c)^2);
end