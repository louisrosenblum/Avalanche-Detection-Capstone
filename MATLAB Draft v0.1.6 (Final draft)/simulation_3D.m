% 3D Simulation
% Louis Rosenblum

%% Initialize program

%clear all
close all

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
    
% Plot signals received by sensors
figure()
hold on
plot(t,signal0); 
plot(t,signal1);
plot(t,signal2);
plot(t,signal3);

plot(t,signal4); 
plot(t,signal5);
plot(t,signal6);
plot(t,signal7);
legend('Sensor 0', 'Sensor 1', 'Sensor 2', 'Sensor 3','Sensor 4', 'Sensor 5', 'Sensor 6', 'Sensor 7');
title("Signals seen by sensors");
xlabel("Time (s)");
ylabel("Amplitude");

%% Run signal processing algorithim

% All possible sensor pair permutations

count = nPr(8,2);

d = [0 1 2 3 4 5 6 7];

p = permn(d,2);

heatmap_out = cell(150,150);

for i = 1:150
    for k = 1:150
        heatmap_out{i,k} = 1;
    end
end

for i = 1:count
    
    if(p(i,1) == p(i,2))
    
    else
    if(p(i,1) == 0)
        ref_sensor = s0;
        ref_signal = signal0;
    elseif(p(i,1) == 1)
        ref_sensor = s1;
        ref_signal = signal1;
    elseif(p(i,1) == 2)
        ref_sensor = s2;
        ref_signal = signal2;
    elseif(p(i,1) == 3)
        ref_sensor = s3;
        ref_signal = signal3;
    elseif(p(i,1) == 4)
        ref_sensor = s4;
        ref_signal = signal4;
    elseif(p(i,1) == 5)
        ref_sensor = s5;
        ref_signal = signal5;
    elseif(p(i,1) == 6)
        ref_sensor = s6;
        ref_signal = signal6;
    else
        ref_sensor = s7;
        ref_signal = signal7;
    end
    
    if(p(i,2) == 0)
        match_sensor = s0;
        match_signal = signal0;
    elseif(p(i,2) == 1)
        match_sensor = s1;
        match_signal = signal1;
    elseif(p(i,2) == 2)
        match_sensor = s2;
        match_signal = signal2;
    elseif(p(i,2) == 3)
        match_sensor = s3;
        match_signal = signal3;
    elseif(p(i,2) == 4)
        match_sensor = s4;
        match_signal = signal4;
    elseif(p(i,2) == 5)
        match_sensor = s5;
        match_signal = signal5;
    elseif(p(i,2) == 6)
        match_sensor = s6;
        match_signal = signal6;
    else
        match_sensor = s7;
        match_signal = signal7;
    end
    
        [heatmap_out] = predict(heatmap_out,f,ref_signal,match_signal,grid,ref_sensor,match_sensor,speed_of_sound);
        
    end
end

heatmap_final = heatmap_out;

%% Draw heatmap
x = [];
y = [];
z = [];

w = [];
compare = 0;

all = [];

for i = 1:150
    for k = 1:150
        x = [x grid{i,k}(1)];
        y = [y grid{i,k}(2)];
        j = heatmap_final{i,k};
        
        all = [all j];
        
        if(j > compare)
            compare = j;
            w = [i k];
        end
        z = [z j];
        
    end
end

figure()
scatter(x,y,1,z)
p = colorbar();
ylabel(p,'Amplitude')

title("Confidence Engine Prediction");
xlabel("X Location (m)");
ylabel("Y Location (m)");

%% Create contour plot

x = 0:1000/328:1000;
y = 0:1000/164:1000;

figure()
[C,h] = contour(x,y,elevation,16);

% Label countours with elevations
%clabel(C,h)
colormap default

hold on


p0 = gscatter(s0(1),s0(2),'Sensor 0'),text(s0(1)+15,s0(2),int2str(s0(3)));
gscatter(s1(1),s1(2),'Sensor 1'),text(s1(1)+15,s1(2),int2str(s1(3)));
gscatter(s2(1),s2(2),'Sensor 2'),text(s2(1)+15,s2(2),int2str(s2(3)));
gscatter(s3(1),s3(2),'Sensor 3'),text(s3(1)+15,s3(2),int2str(s3(3)));
gscatter(s4(1),s4(2),'Sensor 4'),text(s4(1)+15,s4(2),int2str(s4(3)));
gscatter(s5(1),s5(2),'Sensor 5'),text(s5(1)+15,s5(2),int2str(s5(3)));
gscatter(s6(1),s6(2),'Sensor 6'),text(s6(1)+15,s6(2),int2str(s6(3)));
gscatter(s7(1),s7(2),'Sensor 7'),text(s7(1)+15,s7(2),int2str(s7(3)));

p1 = scatter([origin(1)],[origin(2)],'filled','g');

x = grid{w(1),w(2)};
p2 = scatter(x(1),x(2),'filled','b');

xlim([-100 1100]),ylim([-1200 1100]);
title("Spatial Layout");
xlabel("X Position (meters)");
ylabel("Y Position (meters)");

%% Error display

% Plot 100 m error circle

pos = [(origin(1)-100) (origin(2)-100) 200 200]; 
rectangle('Position',pos,'Curvature',[1 1])


avg = mean(all);
dev = std(all);

max = compare;

max_z = (max - avg)/dev;

if(max_z >= 5.894)
    prob = normcdf(max_z,4.5797,0.385) * 100;
    fprintf('Confidence threshold exceeded, system is %d',prob)
    disp('% confident infrasonic signal present')
else
    prob = (1 - (normcdf(max_z,75.57,50.75)))*100;
    fprintf('Confidence below threshold, system is %d',prob)
    disp('% confident infrasonic signal absent');
end

legend([h p0 p1 p2], 'Elevation','Sensor Array','Actual Origin','Ultimate Algorithim Prediction');

hold off

%% Error calculation

prediction = grid{w(1),w(2)};

error = dist2d(origin,prediction)

%% Signal Processing Algorithim Definition

function [heatmap] = predict(heatmap_in,f,signal0,signal1,grid,s0,s1,speed)
    
    x = [];
    y = [];
    z = [];

    
    heatmap = cell(150,150);
   
  % Iterate through all grid points
    for i = 1:150
        for k = 1:150
            
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
        

            signal1_shift = circshift(signal1,round(-shift_1.*f/10));

            % Sum all four signals
            beamformed = signal0 .*signal1_shift;
            
            % Calculate root mean square ampltitude
            amplitude = mean(sqrt(beamformed.^2));
            %amplitude = max(beamformed);
            
            %x = [x grid{i,k}(1)];
            %y = [y grid{i,k}(2)];
            %z = [z amplitude];
            

            
            
            heatmap{i,k} = amplitude .* heatmap_in{i,k};
                        
                        
        end
    end
    
    %figure()
    %scatter(x,y,1,z)
    %p = colorbar();
    %ylabel(p,'Amplitude')

    %title("Confidence Engine Prediction");
    %xlabel("X Location (m)");
    %ylabel("Y Location (m)");
    

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

%% Number of possible permutations function definition

function g = nPr(n,r)

    g = nchoosek(n,r)*factorial(r);
    
end

%% Possible permutatins function

% https://www.mathworks.com/matlabcentral/fileexchange/7147-permn

function [M, I] = permn(V, N, K)

narginchk(2, 3) ;

if fix(N) ~= N || N < 0 || numel(N) ~= 1
    error('permn:negativeN','Second argument should be a positive integer') ;
end
nV = numel(V) ;

if nargin==2 
    
    if nV == 0 || N == 0
        M = zeros(nV, N) ;
        I = zeros(nV, N) ;
    elseif N == 1

        M = V(:) ;
        I = (1:nV).' ;
    else

        [Y{N:-1:1}] = ndgrid(1:nV) ;
        I = reshape(cat(N+1, Y{:}), [], N) ;
        M = V(I) ;
    end
else

    nK = numel(K) ;
    if nV == 0 || N == 0 || nK == 0
        M = zeros(numel(K), N) ;
        I = zeros(numel(K), N) ;
    elseif nK < 1 || any(K<1) || any(K ~= fix(K))
        error('permn:InvalidIndex','Third argument should contain positive integers.') ;
    else
        V = reshape(V, 1, []) ;
        nV = numel(V) ;
        Npos = nV^N ;
        if any(K > Npos)
            warning('permn:IndexOverflow', ...
                'Values of K exceeding the total number of combinations are saturated.')
            K = min(K, Npos) ;
        end
             

        B = nV.^(1-N:0) ;
        I = ((K(:)-.5) * B) ; 
        I = rem(floor(I), nV) + 1 ;
        M = V(I) ;
    end
end

end



