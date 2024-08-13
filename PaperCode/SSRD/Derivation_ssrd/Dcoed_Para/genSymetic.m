function out = genSymetic(ML)

N=length(ML);
inxx=solveRoots(N);
n0=0;
I=1:N;
Im=zeros(inxx,inxx);
for i=1:inxx
    n=solveRoots(N-n0);
    Im(i,i:inxx)=I(n0+1:n+n0);
    Im(i:inxx,i)=I(n0+1:n+n0).';
    n0=n0+n;
end

out=ML(Im);
end

