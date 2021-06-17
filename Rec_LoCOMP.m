function x=Rec_LoCOMP(dic,coff,N)
%% Code by Liu He （aresmiki@163.com）
%% 单原子的重构
np=length(dic);
x=zeros(N,1);
N=min(N,length(coff));
for i=1:N-np+1
    if coff(i)~=0
        x(i:i+np-1)=x(i:i+np-1)+coff(i)*dic(:);
    end
end
end