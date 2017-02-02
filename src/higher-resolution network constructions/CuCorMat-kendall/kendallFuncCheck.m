N = 25218;
L = 1200;
L2 = L * (L - 1) / 2.0;
batchSize = 512;
% threshold = 0.203679;
ep = 1e-6;

%read preprocess results from VS.
% fid = fopen('E:\project\data\4mm\weighted2\103414_resliced_25218_keandall_afterPreprocess.matrix');
% devBold_preprocess = fread(fid, L2 * batchSize ,'float32');
% dp = (reshape(devBold_preprocess,L2, batchSize));

%calculate preprocess results in Matlab.
mp = zeros(L2, batchSize);
numblock = ceil(N/batchSize);
zero = zeros(L,numblock * batchSize - N);
T = [T,zero]; %ceil
addr = 1 + batchSize * 49;
bulk = T(1:L, addr:addr + batchSize - 1);
%matlab preprocess
for i = 1:batchSize
    cnt = 1;
    for j = 1:L
        if(j<L)
            for k = j+1:L
               mp(cnt,i) = sign( bulk(k,i) - bulk(j,i) );
               cnt = cnt + 1;
            end
        end
    end
end