close all;
clear all;

n = 1000; signalPresentAbsent = rand(1,n);
signalPresentAbsent = round(signalPresentAbsent);

for i = 1:length(signalPresentAbsent)
  % if signal present trial
  if signalPresentAbsent(i) == 1
    % then pull a random draw from the signal distribution with mean = 1 and std = 1
    signal(i) = random('norm',1,1);
  else
    % otherwise it is a noise trial so pull a random draw from the noise distribution with mean = 0 and std = 1
    signal(i) = random('norm',0,1);
  end
end

hist(signal);

hist(signal(signalPresentAbsent==1));
% show signal absent distribution

hist(signal(signalPresentAbsent==0));

response = signal>0.5;

% get total number of present trials
nPresent = sum(signalPresentAbsent==1);
% compute hits as all the responses to trials in which signal was present (signalPresentAbsent==1) in which the response was present (i.e. == 1). Divide by number of present trials.
hits = sum(response(signalPresentAbsent==1)==1)/nPresent
% misses are the same except when the responses are 0 (absent even though signal was present)
misses = sum(response(signalPresentAbsent==1)==0)/nPresent
% same idea for correctRejects and falseAlarms
nAbsent = sum(signalPresentAbsent==1);
correctRejects = sum(response(signalPresentAbsent==0)==0)/nAbsent
falseAlarms = sum(response(signalPresentAbsent==0)==1)/nAbsent

zHits = icdf('norm',hits,0,1)
zFalseAlarms = icdf('norm',falseAlarms,0,1)
dPrime = zHits-zFalseAlarms