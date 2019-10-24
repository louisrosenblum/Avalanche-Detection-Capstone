% Team 20 - Avalanche Detection
% Nov 12th, algorithim demo
% Louis Rosenblum, Cayden Seiler, Khristian Jones
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

s0 = [0 0];
s1 = [100 0];
s2 = [0 100];
s3 = [100 100];

%% Grid design

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

%Distance from origin to sensors
d0 = distance(s0,origin) 
d1 = distance(s1,origin)
d2 = distance(s2,origin)
d3 = distance(s3,origin)

dist1 = distance(s0,grid{90,30})
delta1 = d1 - d0
delta2 = d2 - d0
delta3 = d3 - d0

%% Gaussian noise


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