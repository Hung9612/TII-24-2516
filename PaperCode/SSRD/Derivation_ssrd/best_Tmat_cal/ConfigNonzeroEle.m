function [vect,vece]=ConfigNonzeroEle(Tc,Ev)

    for i=1:length(Tc)
        [vect{i}.val,vect{i}.inx] = NonZerosEle(Tc{i},4);
        [vece{i}.val,vece{i}.inx] = NonZerosEle(Ev{i},4); 
    end
end

