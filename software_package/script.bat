echo unweightednetworks
.\exefiles\CUCorMatHr.exe E:/revised/data/4mm E:/revised/data/4mm/mask/mask.nii 0.1 n r p 0.5
.\exefiles\CUCorMat.exe E:/revised/data/4mm E:/revised/data/4mm/mask/mask.nii 0.1 yn n r 0.5
.\exefiles\CUBFW_Lp.exe E:/revised/data/4mm/unweighted 2 gl
.\exefiles\Cp.exe E:/revised/data/4mm/unweighted 2 gk
.\exefiles\SmallWorldProperty.exe E:/revised/data/4mm/unweighted 
.\exefiles\Louvain_Modularity.exe E:/revised/data/4mm/unweighted 2
.\exefiles\PC_CPU.exe E:/revised/data/4mm/unweighted 
.\exefiles\Degree.exe E:/revised/data/4mm/unweighted
.\exefiles\CUBC.exe E:/revised/data/4mm/unweighted
.\exefiles\CUEC.exe E:/revised/data/4mm/unweighted
.\exefiles\ConvertNII.exe E:/revised/data/4mm/unweighted E:/revised/data/4mm/mask/mask.nii 0.1
