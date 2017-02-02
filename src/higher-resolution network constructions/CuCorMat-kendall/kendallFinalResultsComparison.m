
% vsBatch = spa1002(1:128,257:384);
% TransKendall(TransKendall<=0.381512) = 0;

%This is a violent comparison method to find variance
%fighting! never walk backwards

threshold = 0.381512;
batchSize = 128;

for i = 14 :batchSize
    for j = i:batchSize
        str = ['Compare the batch located in ' num2str(i) 'th row and ' num2str(j) 'th column..'];
        disp(str); %keep in mind disp only have a single parameter. 
        %phase 1 partition matrix derived from VS.
        upleft = 1 + batchSize * i;
        downleft = batchSize + batchSize * i;
        columnleft = 1 + batchSize * j;
        columnright = batchSize + batchSize * j;
        vsBatch = spa1002(upleft : downleft, columnleft:columnright);

        %phase 2 calculate and threshold matrix derived from MATLAB.
        TBlock1 = T(:,(columnleft:columnright));
        TBlock2 = T(:,(upleft:downleft));
        
        tic;
        TmKendall = corr(TBlock1, TBlock2, 'type', 'Kendall');
        toc;
        
        TmKendall( isnan(TmKendall) ) = 0;
        TmKendall( TmKendall < 0) = 0;
        TmKendall( TmKendall >= 1 + ep) = 0;
        TmKendall( TmKendall <= threshold) = 0;
        if(i==j)
            TmFinalKendall = TmKendall-diag(diag(TmKendall));
        else
            TmFinalKendall = TmKendall'; %lower triangle stored in row-major format needed to be converted into upper  triangle stored in row-major format
        end

        %phase 3 compare and output
        aggregateError = sum(sum(abs(TmFinalKendall-vsBatch)));
        if(aggregateError<1e-5)
            disp('Matched.');
        else
            str = ['Error occurred! Error is ' num2str(aggregateError) '.'];
            disp(str);
            break;
        end
        
    end
    if(aggregateError>1e-6)
        break;
    end
end
