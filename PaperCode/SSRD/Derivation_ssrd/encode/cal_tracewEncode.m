function [tw,dtw] = cal_tracewEncode(Tmat,pn)
%CALTRACEW 
%   ������α任����
n=length(Tmat);
base=funcBase(n);
% tw=cell(1,n);
dtw=cell(1,n);

for p=1:n
    % ����ÿ����α任����ĺ�������
        [rj,pjw]=conformSgn(Tmat{1,p}.W);
        % �Ե�����α任������б���
%         tw{1,p}.tw=encoders(pjw,base,rj,pn);
        tw=encoders(pjw,base,rj,pn);
    for j=1:p
        dtw{j,p}.difftw=cal_derivite(tw,j,pn);% �Ե�p���任���������еĵ���
    end
end

end

