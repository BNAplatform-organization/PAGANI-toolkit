N = 58523;
L = 1024;
%threshold = 0.372829; for L = 120 and N = 500
%threshold = 0.378431;%new
% batchSize = 128;
% ep = 1e-6;

fid = fopen('E:\project\data\4mm\weightedForNewman\103414_resliced.nii_BOLD_t.matrix');
temp = fread(fid,N * L,'float32');
T = (reshape(temp,L, N));

tic;
kendall = corr(T, T, 'type', 'Pearson');%between each pair of columns in the N-by-P matrix X
toc;
% due to slow calculation, only generate first line;
% kendallFirstLine = zeros(1,N);
% tic;
% for i = 1:10000
%     kendallFirstLine(i) = corr(T(:,4), T(:,i), 'type', 'Kendall'); 
% end
% toc;

% Tbatch2 = T(:,(batchSize+1+batchSize:batchSize+batchSize+batchSize));
% diagflag = false;
% %Tbatch2 = T(:,(1:batchSize));
% Tbatch1 = T(:,(1:batchSize));
% tic;
% mKendall = corr(Tbatch2, Tbatch1, 'type', 'Kendall');
% toc;
% mKendall( isnan(mKendall) ) = 0;
% mKendall( mKendall< 0) = 0;
% mKendall( mKendall >= 1 + ep) = 0;
% if(diagflag == true)
%     mKendall = mKendall-diag(diag(mKendall));
% else
%    TransKendall = mKendall'; 
% end