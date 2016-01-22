function [ r,c,n ] = readcsr( filename )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
fid=fopen(filename);
rlength=fread(fid,1,'int32');
r=fread(fid,rlength,'int32');
clength=fread(fid,1,'int32');
c=fread(fid,clength,'int32');
n=rlength-1;
end

