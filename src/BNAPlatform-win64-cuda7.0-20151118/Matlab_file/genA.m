 [ r,c,n ] = readcsr( 'E:\BNAlauncher\test_new\unweighted\N0001_25218_spa1.000%_cor0.504.csr');
   A  = csr2adjmat(n,r,c);
 
   [ec, D]= IPN_centEigenvector(A);
 