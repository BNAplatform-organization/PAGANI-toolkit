clc;clear;

N = 25218;
L = 1200;
threshold = 0.203679;
ep = 1e-6;

fid = fopen('E:\project\data\4mm\weighted2\103414_resliced.nii_BOLD_t.matrix');
temp = fread(fid,N * L,'float32');
T = (reshape(temp,L, N));
tic;
spearman = corr(T, T, 'type', 'Pearson');%between each pair of columns in the N-by-P matrix X
spearman( spearman<threshold ) = 0;
spearman( spearman>1 ) = 0;
spearman(isnan(spearman)) = 0;
finalResults = spearman-diag(diag(spearman));
toc;

%   read data from disk and check out
   addpath('E:\project\BNAPlatform-win64-cuda7.0-20160115\Matlab_file');
filename = 'E:\project\data\4mm\weighted2\103414_resliced_25218_spearman_spa1.000%_cor0.204_weighted.csr';
[r, c, v, n] = readcsrWeighted(filename);
  A  = csr2adjmatWeighted(n,r,c,v);

  for i=1:N
    for j=1:N
    if(abs(finalResults(i,j)-A(i,j))>ep&& roundn(abs(finalResults(i,j)-A(i,j)),-4) ~= 0.2037  )% error from thrshold;error from storage or computation? 
        disp(i)
        disp(j)
        str=['Matlab results is ' num2str(finalResults(i,j)) '.'];   
        disp(str)
        str=['VS2012 output results is ' num2str(A(i,j)) '.'];   
        disp(str)
        break;
    end
    end
    if(abs(finalResults(i,j)-A(i,j))>ep)
        break;
    end
  end

