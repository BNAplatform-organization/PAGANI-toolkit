addpath('E:\project\BNAPlatform-win64-cuda7.0-20160115\Matlab_file');
filename = 'E:\project\data\4mm\weighted2\103414_resliced_25218_kendall_spa42.504%_cor0.000_weighted.csr';
[r, c, v, n] = readcsrWeighted(filename);
 Ashield  = csr2adjmatWeighted(n,r,c,v);
 
 % read fcs
 