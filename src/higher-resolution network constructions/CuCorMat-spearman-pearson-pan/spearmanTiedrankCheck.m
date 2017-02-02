L = 1200;
N = 25218;

Rank = zeros(L, N);
tic;
for i = 1:N
Rank(:,i) = tiedrank(T(:,i));
end
toc;

%load rankfile from vs2012
fid = fopen('E:\project\data\4mm\weighted2\103414_resliced.nii_BOLD_t_AfterRank.matrix');
tempAfterRank = fread(fid,N * L,'float32');
RankVS = (reshape(tempAfterRank,L, N));
difference = RankVS - Rank;

%verify the formula

single = T(:,4);
singleRank = zeros(1,L);
for i = 1:L
    benchmark = single(i);
    lessOne = 0;
    equalOne = 0;
for j = 1:L
    temp = single(j);
    if(temp < benchmark)
        lessOne = lessOne + 1;
    end
    if( temp == benchmark )
        equalOne = equalOne + 1 ;
    end
end
    singleRank(i) = lessOne + 0.5 * ( 1 + equalOne );
end

