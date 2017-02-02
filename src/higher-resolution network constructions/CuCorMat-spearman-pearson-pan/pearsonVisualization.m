clear;
fid = fopen('E:\project\data\4mm\weightedPearson\0_5\103414_resliced_25218_pearson_spa0.500%_cor0.240_weighted_deg.nm');
length1 = fread(fid,1,'int32');
deg1 = fread(fid,length1,'float32');

fid = fopen('E:\project\data\4mm\weightedPearson\0_5\105115_resliced_25218_pearson_spa0.500%_cor0.244_weighted_deg.nm');
length2 = fread(fid,1,'int32');
deg2 = fread(fid,length2,'float32');

fid = fopen('E:\project\data\4mm\weightedPearson\0_5\115320_resliced_25218_pearson_spa0.500%_cor0.267_weighted_deg.nm');
length3 = fread(fid,1,'int32');
deg3 = fread(fid,length3,'float32');

fid = fopen('E:\project\data\4mm\weightedPearson\0_5\117122_resliced_25218_pearson_spa0.500%_cor0.301_weighted_deg.nm');
length4 = fread(fid,1,'int32');
deg4 = fread(fid,length4,'float32');

if(~(length1==length2 && length2==length3 && length3==length4))
    disp('error!')
else
    deg5 = zeros(length1,1);
    for i=1:length1
        temp = [deg1(i,1) deg2(i,1) deg3(i,1) deg4(i,1)];
        deg5(i,1) = mean(temp);
    end
end

fid = fopen('E:\project\data\4mm\weightedPearson\average\25218_pearson_spa0.500%_deg_aver.nm','w');
fwrite(fid, length1, 'int32')
fwrite(fid, deg5, 'float32')
fclose(fid);

