%% Team 20 - Avalanche Detection
% Nov 12th, algorithm demo
% Louis Rosenblum, Cayden Seiler, Khristian Jones
<<<<<<< Updated upstream
%% Sensor and Grid locations

=======
%% Noise for Waveform
fileReader = dsp.AudioFileReader('Avy.wav');
writer = audioDeviceWriter('SampleRate', fileReader.SampleRate);
scope = dsp.TimeScope(1,...
                      fileReader.SampleRate, ...
                      'TimeSpanOverrunAction','Scroll', ...
                      'TimeSpan',6.5, ...
                      'BufferLength',1.5e6, ...
                      'YLimits',[-1 1],...
                      'ShowGrid',true,...
                      'ShowLegend',true);
%% Sensor placement
>>>>>>> Stashed changes



%% Generate random temps/humidity/air pressure



%% Calc speed of sound


%% Randomly pick origin


%% Calculate distance to sensor



%% Generate phase delays



<<<<<<< Updated upstream
%% 4 copies of .wav file with delays



%% Add noise 


=======
%% Calculate distance to sensors

<<<<<<< Updated upstream
%Distance from origin to sensors
d0 = distance(s0,origin) 
d1 = distance(s1,origin)
d2 = distance(s2,origin)
d3 = distance(s3,origin)

dist1 = distance(s0,grid{90,30})
delta1 = d1 - d0
delta2 = d2 - d0
delta3 = d3 - d0
>>>>>>> Stashed changes


<<<<<<< Updated upstream
%% The meat and potatoes of this bad boi
%  Calculate origin of signal
=======
d0 = distance(s0,origin);
d1 = distance(s1,origin);
d2 = distance(s2,origin);
d3 = distance(s3,origin);

delta1 = d1 - d0;
delta2 = d2 - d0;
delta3 = d3 - d0;


%% Signal Generation

figure();
t = 0:1/3413:0.3;


signal0 = cos(10*2*pi.*t);

wavelength = speed_of_sound/10;
shift1 = delta1/wavelength;
shift2 = delta2/wavelength;
shift3 = delta3/wavelength;

signal1 = cos(10*2*pi.*(t-shift1/10));
signal2 = cos(10*2*pi.*(t-shift2/10));
signal3 = cos(10*2*pi.*(t-shift3/10));

% Add gaussian noise

signal0 = awgn(signal0,25);
signal1 = awgn(signal1,25);
signal2 = awgn(signal2,25);
signal3 = awgn(signal3,25);

plot(t,signal0), hold on
plot(t,signal1);
plot(t,signal2);
plot(t,signal3);
% legend('Sensor 0', 'Sensor 1', 'Sensor 2', 'Sensor 3');
% title("Signals seen by sensors");
% xlabel("Time");
% ylabel("Amplitude"); %hold off;

<<<<<<< Updated upstream
% filt0 = lowpass(signal0,20,10000);
% filt1 = lowpass(signal1,20,10000);
% filt2 = lowpass(signal2,20,10000);
% filt3 = lowpass(signal3,20,10000);
% 
% plot(t,filt0), hold on
% plot(t,filt1);
% plot(t,filt2);
% plot(t,filt3);
% legend('Sensor 0', 'Sensor 1', 'Sensor 2', 'Sensor 3');
% title("Signals seen by sensors");
% xlabel("Time");
% ylabel("Amplitude"); hold off;
=======
filt0 = lowpass(signal0,20,10000);
filt1 = lowpass(signal1,20,10000);
filt2 = lowpass(signal2,20,10000);
filt3 = lowpass(signal3,20,10000);

plot(t,filt0,'LineWidth',1.5), %hold on
plot(t,filt1,'LineWidth',1.5);
plot(t,filt2,'LineWidth',1.5);
plot(t,filt3,'LineWidth',1.5);
legend('Sensor 0', 'Sensor 1', 'Sensor 2', 'Sensor 3','Filt 0','Filt 1','Filt 2','Filt 3');
title("Filtered Signals");
xlabel("Time");
ylabel("Amplitude"); hold off;

>>>>>>> Stashed changes

amplitude = max(signal0(:));


%% Algorithim execution

<<<<<<< Updated upstream
%[guess, height] = algorithm(s0,s1,s2,s3,signal0,signal1,signal2,signal3,grid,speed_of_sound,0);
[filteredGuess, height2] = algorithm(s0,s1,s2,s3,signal0,signal1,signal2,signal3,grid,speed_of_sound);
=======
[guess, height] = algorithm(s0,s1,s2,s3,filt0,filt1,filt2,filt3,grid,speed_of_sound);

>>>>>>> Stashed changes

%% Plot 

figure();
% Sensors
%gscatter([0 100 0 100],[0 0 100 100],[0;1;2;3]),
gscatter(0,0,'Sensor 0', 'b'),hold on
gscatter(0,100,'Sensor 1', 'r');
gscatter(100,0,'Sensor 2', 'y');
gscatter(100,100,'Sensor 3', 'm');
xlim([-100 1100]),ylim([-100 2100]);


% True origin
scatter([origin(1)],[origin(2)],'filled');
%scatter([guess(1)],[guess(2)],'filled');
scatter([filteredGuess(1)],[filteredGuess(2)],'filled');
legend('Sensor 0', 'Sensor 1', 'Sensor 2', 'Sensor 3', 'Origin','Predicted Origin','Filtered Predicted Origin');
title("Sensor Grid"); 
%xlabel("X (m)");
%ylabel("Y (m)");
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
    plot(x1,y1,'b','HandleVisibility', 'off'), hold on;
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
    plot(x1,y1,'b','HandleVisibility','off'),xlabel("m"),ylabel("m")
    x1 = [];
    y1 = [];
end
hold off;

%% Error calculation

d_1 = distance(s0,origin);
d_2 = distance(s0,guess);

percent_error = sqrt((d_2 - d_1)^2)/d_1 * 100

%% Prediction algorithm

function [predict, amp] = algorithm(s0,s1,s2,s3,signal_0,signal_1,signal_2,signal_3,grid,speed)
    
    amp = 0;
    predict = {1,1};
    for i = 1:100
        for k = 1:100
            distance0 = distance(s0,grid{i,k});
            distance1 = distance(s1,grid{i,k});
            distance2 = distance(s2,grid{i,k});
            distance3 = distance(s3,grid{i,k});
            
            delta_1 = distance1 - distance0;
            delta_2 = distance2 - distance0;
            delta_3 = distance3 - distance0;
            
            wave_length = speed/10;
            
            shift_1 = delta_1/wave_length;
            shift_2 = delta_2/wave_length;
            shift_3 = delta_3/wave_length;
            
            signal1_shift = circshift(signal_1,round(-shift_1*1024/3));
            signal2_shift = circshift(signal_2,round(-shift_2*1024/3));
            signal3_shift = circshift(signal_3,round(-shift_3*1024/3));
            
            beamformed = signal_0 + signal1_shift + signal2_shift + signal3_shift;
            
            % Root mean square ampltitude

                filtered = lowpass(beamformed,20,10000);
                filtered = (filtered).^2;
                filtered = sqrt(filtered);
                amplitude = mean(filtered);
            
                 %beamformed = (beamformed).^2;
                 %beamformed = sqrt(beamformed);
                 %amplitude = mean(beamformed);
            
            if amplitude > amp
                amp = amplitude;
                predict = grid{i,k};
            end
            
        end
    end
end
>>>>>>> Stashed changes



%% Plot sensors, grid, predicted origin, actual origin




%% Print if matching. If not how far off is it.
=======

while ~isDone(fileReader)
    wave1 = fileReader();
    %wave2 = fileReader();
    %wave3 = fileReader();
    %wave4 = fileReader();
    wc1 = wave1 + (2e-1/4) * randn(1024,1);
    %wc2 = wave2 + (1e-2/4) * randn(1024,1);
    %wc3 = wave3 + (1e-2/4) * randn(1024,1);
    %wc4 = wave4 + (1e-2/4) * randn(1024,1);
    writer(wc1)
    scope([wc1, wave1]);
end
release(fileReader)
release(scope)
release(writer)

%% Distance function definition

function dist = distance(p1,p2)
    a = p2(1);
    b = p2(2);
    dist = sqrt(abs((p2(1) - p1(1))^2 + (p2(2)-p1(2))^2));
end
>>>>>>> Stashed changes
