addpath('E:\project\BNAPlatform-win64-cuda7.0-20160115\Matlab_file');
filename = 'E:\project\data\4mm\weighted2\103414_resliced_25218_kendall_spa1.002%_cor0.382_weighted.csr';
[r, c, v, n] = readcsrWeighted(filename);
spa1002  = csr2adjmatWeighted(n,r,c,v);

% threshold = 0.378431;
% A = A - diag(diag(A));
% sum(sum(A~=0))