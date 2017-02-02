sampleSize = 128;
dimension = 1200000;
practical = dimension + 200000;

%x = rand(sampleSize, dimension);
x = rand(practical,sampleSize);

fid = fopen('E:\project\CuCorMat-spearman-pearson-pan\Data.datt','w');
fwrite(fid,x,'float32');
fclose(fid);