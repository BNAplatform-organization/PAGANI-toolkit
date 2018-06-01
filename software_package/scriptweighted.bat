echo weightednetworks
.\exefiles_weighted\CUCorMatHr.exe E:/data/testdata E:/data/testdata/mask.nii 0.2 n s p 0.1
.\exefiles_weighted\Cp.exe E:/data/testdata/weighted 0 2  gn
.\exefiles_weighted\Degree.exe E:/data/testdata/weighted
.\exefiles_weighted\CUEC.exe E:/data/testdata/weighted
.\exefiles_weighted\ConvertNII.exe E:/data/testdata/weighted E:/data/testdata/mask.nii 0.2
