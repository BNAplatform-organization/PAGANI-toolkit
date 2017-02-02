N = 25218; 
for i=1:N
    for j=1:N
    if(abs(Aerror(i,j)-YuanShiJuZhen(i,j))>ep) %&& Aerror(i,j)~=0 )
        disp(i)
        disp(j)
        str=['Matlab results is ' num2str(YuanShiJuZhen(i,j)) '.'];   
        disp(str)
        str=['VS2012 output results is ' num2str(Aerror(i,j)) '.'];   
        disp(str)
%         break;
    end
    end
%     if(abs(Aerror(i,j)-A(i,j))>ep) %&& Aerror(i,j)~=0 )
%         break;
%    end
  end