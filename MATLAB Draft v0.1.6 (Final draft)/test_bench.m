diary('log.txt')

result = []
prob = []

z_scores = []
probs = []
tic



for snr = 8:30
for go = 1:10
    snr
    simulation_3D();
    result(go) = max_z;
    prob(go) = prob;
end
    z_scores(snr-7) = mean(result);
    probs(snr - 7) = mean(prob);
end


diary('off')

toc

