echo weightednetworks
.\exefiles_weighted\CUCorMatHr.exe E:/revised/data/4mm E:/revised/data/4mm/mask.nii 0.1 n s p 0.5 
.\exefiles_weighted\CUBFW_Lp.exe E:/revised/data/4mm/weighted 2 gnl
.\exefiles_weighted\Cp.exe E:/revised/data/4mm/weighted 2  1  gnk
.\exefiles_weighted\SmallWorldProperty.exe E:/revised/data/4mm/weighted 
.\exefiles_weighted\Louvain_Modularity.exe E:/revised/data/4mm/weighted 2
.\exefiles_weighted\PC_CPU.exe E:/revised/data/4mm/weighted 
.\exefiles_weighted\Degree.exe E:/revised/data/4mm/weighted
.\exefiles_weighted\CUEC.exe E:/revised/data/4mm/weighted
.\exefiles_weighted\ConvertNII.exe E:/revised/data/4mm/weighted E:/revised/data/4mm/mask.nii 0.1
