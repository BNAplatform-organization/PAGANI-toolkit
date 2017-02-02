batchSize = 6144;
ep = 1e-6;

% fid = fopen('E:\project\data\4mm\weighted2\103414_resliced_25218_pearsoncompareR_right.matrix');%true
% devCormatTrue = fread(fid, batchSize * batchSize ,'float32');
% dcTrue = (reshape(devCormatTrue,batchSize, batchSize));%reshape don,t change representation form for column-major storage case.

fid = fopen('E:\project\data\1mm\103414_1mm.nii\weighted2\103414_1mm_1561152_spearmancompareW_wrong.matrix');%false
devCormatRight = fread(fid, batchSize * batchSize ,'float32');
dcWrong = (reshape(devCormatRight,batchSize, batchSize));


% for row = 1: batchSize
%     for column = 1:batchSize
%         if( abs(dc(row, column) - 0.628851532936096)<ep )
%             disp(row)
%             disp(column)
%             break;
%         end
%     end
%      if( (dc(row, column) - 0.628851532936096)<ep )
%         break;
%      end
% end
