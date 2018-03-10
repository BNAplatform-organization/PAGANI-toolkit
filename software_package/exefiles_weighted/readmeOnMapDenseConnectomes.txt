This is the specifications regarding network construction. Three kinds of functional connectivity measures were intergreted in PAGANI toolkit, including Pearson, Kendall, Spearman coefficient of correlation.
Several scripts are provided as name suggests. It's recommendable that users should call network construction functions by altering and running these scripts. 
Generally, a standard command template is:
.\CuCorMat-sca-kendall.exe   E:\project\data\4mm E:\project\data\BNU_data_S2_proto-4mm\mask.nii 0.2 n s p 0.1 1 10
respcetively to correspond meaning as follows:
procedure_name data_directory  mask_file_path  mask_threshold  to_average  network_threshold_type coefficient_type network_threshold_value                                 
Interpretation:
1. procedure_name implies the type of functional connectivity measures you want established newtorks to have. For CuCorMat-sca-Pspear, Pearson and Spearman; For CuCorMat-sca-kendall, Kendall.
1. data_directory is a directory where all files with .nii format will be tackled as the input files of the program.
2. mask_file_path is the path of mask file. Laying mask file in data_directory is generally not recommended since error(s) may occur(s) from time to time.
3. mask_threshold: You guys konw it as name suggests!
4. to_average: this parameter can only be "n" currently. Please dont reason why..
5. network_threshold_type: alternative with "r" or "s". "r" for correlation-value-based thresholding method, and "s" for sparsity-based thresholding method.
6. coefficient_type: "s" for only calculating spearman coefficient of correlation; "p" for for only calculating pearson coefficient of correlation; "sp" or "ps" for calculating both correlation coefficient. Of note,for procedure 'CuCorMat-sca-kendall', no this parameter.
7. network_threshold_value: predicated on parameter network_threshold_type. If network_threshold_type is "r", you can put into multiple values between 0 to 1 spearated by spaces; Else for "s" you can put into multiple values between 0 to 100 without adding sign "%". 
Best regards.