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
s2 = [0 100];
s3 = [100 100];

%% Signal condition generation

% Two random intergers from 1-329 for origin out of possible grid indexes
randx = randi(329,1,1)
randy = randi(329,1,1)

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
snr = randi([16 32])

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

%% Run signal processing algorithim

[w heatmap avg std] = predict(signal0,signal1,signal2,signal3,grid,s0,s1,s2,s3,speed_of_sound);

%% Draw heatmap
x = [];
y = [];
z = [];

for i = 1:329
    for k = 1:329
        x = [x grid{i,k}(1)];
        y = [y grid{i,k}(2)];
        k = (heatmap{i,k});
        z = [z k];
        
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

gscatter(0,0,'Sensor 0', 'b'),text(15,0,"2885 m");
gscatter(0,100,'Sensor 1', 'r'),text(115,0,"2772 m");
gscatter(100,0,'Sensor 2', 'y'),text(15,100,"2648 m");
gscatter(100,100,'Sensor 3', 'm'),text(115,100,"2560 m");

scatter([origin(1)],[origin(2)],'filled');

x = grid{w(1),w(2)};
scatter(x(1),x(2),'filled');

xlim([-100 1100]),ylim([-100 3100]);
title("Spatial Layout");
xlabel("X Position (meters)");
ylabel("Y Position (meters)");

%% Error display

% Plot error circle

pos = [(origin(1)-100) (origin(2)-100) 200 200]; 
rectangle('Position',pos,'Curvature',[1 1])

hold off

legend('Sensor 0', 'Sensor 1', 'Sensor 2', 'Sensor 3','Actual Origin','Highest amplitude alignment');
%% Error calculation

    check = 0;
    error = [];

for i = 1:length(pass)/2
    
    x = pass(2*i-1);
    y = pass(2*i);
    
    point = grid{x,y};
    point = point(1:2);
    
    origin_2d = origin(1:2);
    
    dist_2d = dist2d(point,origin_2d)
    
    if(dist_2d > 100)
        check = 1;
        error = [error (dist_2d^2/(100^2) - 1)*100];
    else
        error = [error 0];
    end
end

    if(check == 1)
       disp("Failure: One or more threshold exceeding predictions over 100m off true origin");
       percent_error = mean(error)
    else
        disp("Success: All threshold exceeding predictions within 100m of true origin");
        percent_error = 0
    end



%% Signal Processing Algorithim Definition

function [a heatmap avg std_mag] = predict(signal0,signal1,signal2,signal3,grid,s0,s1,s2,s3,speed)

    data = [];
    amp = 0;
    
    heatmap = cell(329,329);
   
  % Iterate through all grid points
    for i = 1:329
        for k = 1:329
            
            % Calculate distance from current grid point to each sensor
            distance0 = distance(s0,grid{i,k});
            distance1 = distance(s1,grid{i,k});
            distance2 = distance(s2,grid{i,k});
            distance3 = distance(s3,grid{i,k});
            
            % Determine difference in distance to reach sensor 1-3 compared
            % to reference sensor 0
            delta_1 = distance1 - distance0;
            delta_2 = distance2 - distance0;
            delta_3 = distance3 - distance0;
            
            % Calculate wavelength from speed of sound
            wave_length = speed/10;
            
            % Calculate phase shifts from wavelength
            shift_1 = delta_1/wave_length;
            shift_2 = delta_2/wave_length;
            shift_3 = delta_3/wave_length;
            
            % Shift signals 1-3 accordingly, in attempt to match signal 0
            
            val1 = round(-shift_1.*4096/12);
            val2 = round(-shift_2.*4096/12);
            val3 = round(-shift_3.*4096/12);
            
            
            signal1_shift = circshift(signal1,round(-shift_1.*4096/12));
            signal2_shift = circshift(signal2,round(-shift_2.*4096/12));
            signal3_shift = circshift(signal3,round(-shift_3.*4096/12));
            
            
  
            
            % Sum all four signals
            beam1 = signal0 .*signal1_shift;
            beam2 = signal0 .*signal2_shift;
            beam3 = signal0 .*signal3_shift;
            
            beamformed = beam1 + beam2 + beam3;
            
            % Calculate root mean square ampltitude
            amplitude = mean(sqrt(beamformed.^2));
            %amplitude = max(beamformed);
            
            if amplitude > amp
                amp = amplitude;
                a = [i k];
                sig0 = signal0;
                sig1 = signal1_shift;
                sig2 = signal2_shift;
                sig3 = signal3_shift;
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
    plot(t,sig2)
    plot(t,sig3)
    
    legend('Sensor 0', 'Sensor 1', 'Sensor 2', 'Sensor 3');
    title("Maximum Peak Alignment");
    xlabel("Time (s)");
    ylabel("Amplitude");

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