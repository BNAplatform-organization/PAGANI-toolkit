function [ A ] = csr2adjmat(n,r,c)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
A=zeros(n);
r = r+1;
c = c+1;
for i=1:n
    for j=r(i):(r(i+1)-1)
        A(i,c(j)) = 1;
    end    
end

end

