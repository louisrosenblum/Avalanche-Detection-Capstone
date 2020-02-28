microphone = phased.CustomMicrophoneElement('FrequencyVector',[20,20000],'FrequencyResponse',[1,1]);
array = phased.ULA('Element',microphone,'NumElements',6,'ElementSpacing',2);
fs = 8000;
t = 0:1/fs:0.3;
x = chirp(t,0,1,500);
c = 340;
collector = phased.WidebandCollector('Sensor',array,...
    'PropagationSpeed',c,'SampleRate',fs,'ModulatedInput',false);
incidentAngle = [-50;30];
x = collector(x.',incidentAngle);

sigma = 0.2;
noise = sigma*randn(size(x));
rx = x + noise;

beamformer = phased.TimeDelayBeamformer('SensorArray',array,...
    'SampleRate',fs,'PropagationSpeed',c,...
    'Direction',incidentAngle);
y = beamformer(rx);

plot(t,rx(:,1),t,y)
xlabel('Time (sec)')
ylabel('Amplitude')
legend('Original','Beamformed')