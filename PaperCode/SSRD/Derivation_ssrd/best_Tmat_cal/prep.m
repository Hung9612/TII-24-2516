function [vect,vece] = prep(Tc,Ev)

n=length(Tc);

for i=1:n
    [vect{i}.val,vect{i}.inx] = NonZerosEle(Tc{i},4);
    [vece{i}.val,vece{i}.inx] = NonZerosEle(Ev{i},4); 
end
end

