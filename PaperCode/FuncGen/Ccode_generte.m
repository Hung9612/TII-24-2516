function  Ccode_generte(rbt)

Theta=rand(length(rbt.Theta),1);
Ptensor = cell2tensor(rbt);
IDMname=sprintf('TorqueCompute_%s.c',rbt.name);

addpath(fullfile(makeDir, 'genfile'));
filename_inverse = fullfile(makeDir, 'genfile', IDMname);
doMex = true; 
GenInverseDynamics(filename_inverse, uint8(rbt.tauw), Ptensor, Theta, 'mex', doMex);
end

function path = makeDir
    path = fullfile( fileparts(mfilename('fullpath')), '..');
end

function Ptensor = cell2tensor(rbt)
    Ptensor=zeros(size(rbt.tauw,2),size(rbt.Theta,1),rbt.DoF);
    for i=1:rbt.DoF
        Ptensor(:,:,i)=round(rbt.Pmap{i},6);
    end
    
end
