function random_perturbation_generator(n, dirlabel, M);

% generates M nxn matrices of pseudo-random entries, uniformly distributed in
% [-1,1], saves them as 'RPnM_dirlabel/RPj.mat' where j=1...M

dirname = sprintf('RP%gn%gM_%s', n, M, dirlabel);
unix(sprintf('mkdir %s', dirname));

for j=1:M,
    RP = 2*(rand(n)-.5);
    savefile = sprintf('%s/RP%g', dirname, j);
    save(savefile, 'RP');
end    

