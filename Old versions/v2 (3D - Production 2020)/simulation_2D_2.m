% 3D Simulation
% Louis Rosenblum

%% Initialize program

clear all
close all

%% Populate grid with X,Y coordinates

grid = cell(329,329);

for i = 1:329
    for k = 1:329
        x = 1000*(i-1)/328;
        y = 1000*(k-1)/328+2000;
        grid{i,k} = [x y];

    end

end

%% Sensor placement

% Hardcoding sensor locations at (x,y) coordinates
s0 = [0 0];
s1 = [100 0];

s2 = [900 0];
s3 = [1000 0];

%% Signal condition generation

% Two random intergers from 1-329 for origin out of possible grid indexes
randx = randi(329,1,1)
randy = randi(329,1,1)

origin = grid{randx, randy};

% Generate random temp in celsius, -40 C to 10 C
tempc = -10;

% Calculate speed of sound in m/s
speed_of_sound = 331.3 * sqrt(1 + (tempc / 273.15))


%% Calculate distance to sensors

d0 = distance(origin,s0);
d1 = distance(origin,s1);
d2 = distance(origin,s2);
d3 = distance(origin,s3);

% Calculate difference in distance from sensors 1-3 to reference sensor 0
delta1 = d1 - d0;
delta2 = d2 - d0;
delta3 = d3 - d0;

% Calculate amplitude decay over each distance based on energy distributed
% over surface area of a sphere
decay0 = 100000000/(4*pi*d0^2);
decay1 = 100000000/(4*pi*d1^2);
decay2 = 100000000/(4*pi*d2^2);
decay3 = 100000000/(4*pi*d3^2);
%% Signal Generation

% Time vector
t = -0.6:1/3413:0.6;

wave = rand(1,length(t));

% Generate signal hitting the reference sensor
signal0 = decay0 .* sinc(10*2*pi.*t) .* heaviside(t).*wave;

% Calculate wavelength
wavelength = speed_of_sound/10;

% Shift other signals to match distance travelled to each sensor
shift1 = delta1/wavelength;
shift2 = delta2/wavelength;
shift3 = delta3/wavelength;

% Generate signals received by each sensor
signal1 = decay1 .* sinc(10*2*pi.*(t-shift1/10)).* heaviside(t-shift1/10).*wave;
signal2 = decay2 .* sinc(10*2*pi.*(t-shift2/10)).* heaviside(t-shift2/10).*wave;
signal3 = decay3 .* sinc(10*2*pi.*(t-shift3/10)).* heaviside(t-shift3/10).*wave;

% Power factor
pf = 0;
snr = 30;

% Add independent gaussian noise to each signal
signal0 = awgn(signal0,snr,pf);
signal1 = awgn(signal1,snr,pf);
signal2 = awgn(signal2,snr,pf);
signal3 = awgn(signal3,snr,pf);
    
% Plot signals received by sensors
figure()
hold on
plot(t,signal0); 
plot(t,signal1);
plot(t,signal2);
plot(t,signal3);
legend('Sensor 0', 'Sensor 1', 'Sensor 2', 'Sensor 3');
title("Signals seen by sensors");
xlabel("Time (s)");
ylabel("Amplitude");

%% Run signal processing algorithim for s0, s1

[heatmap avg std] = predict(signal0,signal1,grid,s0,s1,speed_of_sound);

%% Run again for s2,s3

[heatmap2 avg2 std2] = predict(signal2,signal3,grid,s2,s3,speed_of_sound);

%% Draw heatmap
x = [];
y = [];
z = [];

track = 0;
w = []
for i = 1:329
    for k = 1:329
        x = [x grid{i,k}(1)];
        y = [y grid{i,k}(2)];
        j = (heatmap{i,k}*heatmap2{i,k});
        
        if (j > track)
            track = j;
            w = [i k]
        end
        z = [z j];
        
    end
end

figure()
scatter(x,y,1,z)
p = colorbar();
title(p,'Magnitude')

title("Spatial Alignment Amplitude");
xlabel("X Location (m)");
ylabel("Y Location (m)");

%% Create contour plot

x = 0:1000/328:1000;
y = 2000:1000/164:3000;

figure()

hold on

gscatter(0,0,'Sensor 0', 'b');
gscatter(100,0,'Sensor 1', 'r');
gscatter(900,0,'Sensor 2', 'g');
gscatter(1000,0,'Sensor 3', 'y');

scatter([origin(1)],[origin(2)],'filled');

x = grid{w(1),w(2)};
scatter(x(1),x(2),'filled');

xlim([-100 1100]),ylim([-100 3100]);
title("Spatial Layout");
xlabel("X Position (meters)");
ylabel("Y Position (meters)");

legend('Sensor 0', 'Sensor 1','Signal origin','Algorithim prediction');

%% Signal Processing Algorithim Definition

function [heatmap avg std_mag] = predict(signal0,signal1,grid,s0,s1,speed)

    data = [];
    amp = 0;
    
    x = [];
    y = [];
    
    heatmap = cell(329,329);
   
  % Iterate through all grid points
    for i = 1:329
        for k = 1:329
            
            % Calculate distance from current grid point to each sensor
            distance0 = distance(s0,grid{i,k});
            distance1 = distance(s1,grid{i,k});
            
            % Determine difference in distance to reach sensor 1-3 compared
            % to reference sensor 0
            delta_1 = distance1 - distance0;
            
            % Calculate wavelength from speed of sound
            wave_length = speed/10;
            
            % Calculate phase shifts from wavelength
            shift_1 = delta_1/wave_length;

            % Shift signals 1-3 accordingly, in attempt to match signal 0
            
            val1 = round(-shift_1.*4096/12);
            
            
            signal1_shift = circshift(signal1,round(-shift_1.*4096/12));
            
            % Sum all four signals
            beam1 = signal0 .*signal1_shift;
            
            beamformed = beam1;
            
            % Calculate root mean square ampltitude
            amplitude = mean(sqrt(beamformed.^2));
            %amplitude = max(beamformed);
            
            x = [x shift_1];
            y = [y amplitude];
            
            if amplitude > amp
                amp = amplitude;
                sig0 = signal0;
                sig1 = signal1_shift;
            end
            
            heatmap{i,k} = amplitude;
   
            data = [data amplitude];
            
        end
    end
   
    avg = mean(data);
    std_mag = std(data);
    
    % Time vector
    t = -0.6:1/3413:0.6;
    figure();
    
    plot(t,sig0), hold on
    plot(t,sig1)
    
    legend('Sensor 0', 'Sensor 1');
    title("Maximum Peak Alignment");
    xlabel("Time (s)");
    ylabel("Amplitude");
    
    figure()
    scatter(x,y)
    xlabel('Phase delay')
    ylabel('Aligned amplitude')

end

    
%% Distance function definition

function dist = distance(p1,p2)
    a = p1(1);
    b = p1(2);
    
    d = p2(1);
    e = p2(2);

    dist = sqrt((d-a)^2+(e-b)^2);
end

%% 2D Distance function definition

function dist = dist2d(p1,p2)
    a = p1(1);
    b = p1(2);
    
    d = p2(1);
    e = p2(2);

    dist = sqrt((d-a)^2+(e-b)^2);
end

%% LogBASE function defintion

function k = logbase(b,x)
    k = log(x)/log(b);
end