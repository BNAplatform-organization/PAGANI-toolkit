addpath('E:\project\CuCorMat-spearman-pearson-pan');
sampleSize = 128;
dimension = 58523;
x = rand(sampleSize, dimension);
tic;
r = IPN_fastCorr(x,x);
toc;