clear;
addpath('E:\project\CuCorMat-spearman-pearson-pan');
sampleSize = 1024;
dimension = 58523;

x = rand(sampleSize, dimension);
tic;
R = corr(x, x, 'type', 'Spearman');
toc;