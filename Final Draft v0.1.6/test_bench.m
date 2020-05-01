diary('log.txt')


tic

snr_error = []

for snr = 0.5:0.5:7
    results = [];
for go = 1:10
    snr
    simulation_3D();
    results(go) = error;
end
    snr_error = [snr_error mean(results)];
end

diary('off')

toc

