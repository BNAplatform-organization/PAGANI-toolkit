%note: this file can only work after running kendallFuncCheck.m

batchSize = 1024;
fid = fopen('E:\project\data\4mm\weighted2\103414_resliced_25218_keandall_afterGemm.matrix');%true
devCormat = fread(fid, batchSize * batchSize ,'float32');
dc = (reshape(devCormat,batchSize, batchSize));

%calculate batch in Matlab
%mc = mp' * mp;

fid = fopen('E:\project\data\4mm\weighted2\103414_resliced_25218_keandall_afterPreprocess.matrix'); %false
devCormat_postpreprocess = fread(fid, batchSize * batchSize ,'float32');
dcpp = (reshape(devCormat_postpreprocess,batchSize, batchSize));