function [rateW,Wi] = conformSgn(Wi)
%UNTITLED2 exchange the sign of each function
%   INPUT:
%       -Wi: function, class: syms
%
%   OUTPUT:
%       -Wi:positive function
%       -rateW: sing vector

rateW=ones(1,length(Wi));
for ki=1:length(Wi)
   strWi=char(Wi(ki));
   if ismember('-',strWi(1))
       rateW(ki)= -1;
       Wi(ki)=-1*Wi(ki);
   end
   if Wi(ki)==1||Wi(ki)==0
       % since diff the Transformation matrix w.r.t joint, constant to be
       % zeros
      rateW(ki)=0;
      Wi(ki)=0;
   end
end

end

