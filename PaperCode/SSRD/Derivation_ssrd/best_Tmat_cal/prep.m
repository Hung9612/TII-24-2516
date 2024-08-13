function [vect,vece] = prep(Tc,Ev)
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明
n=length(Tc);

for i=1:n
    [vect{i}.val,vect{i}.inx] = NonZerosEle(Tc{i},4);
    [vece{i}.val,vece{i}.inx] = NonZerosEle(Ev{i},4); 
end
end

