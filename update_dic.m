function d=update_dic(y,x,dic)
%% Code by Liu He （aresmiki@163.com）
% x  是系数， y 是需要近似的函数 ， np 是字典原子大小
%% 用 x构造卷积字典
np=length(dic);
N=size(x(:),1);
%%
% XDic=circulant(x(:));
% XDic=XDic(:,1:np);
XDic=get_Circshift(repmat(x(:),1,np),[0:1:np-1]); 
y=y(:);
%  XDic*dic 是卷积结果. %这部分可直接用卷积代替矩阵乘法
%% argmin ||y-XDic*dic||_F 
% %% 用梯度下降求解
% options = optimoptions(@fmincon,'HessianApproximation','lbfgs','MaxIterations',20);
% d=fmincon(@(xv)Obj_min(y(1:N),XDic,xv),dic,[],[],[],[],[],[],[],options);
% d=cgls(XDic,y(1:N),0,1e-3,20,'true',dic);
d=cgls(XDic,y(1:N),0,1e-3,20,[],dic); %共轭梯度下降
% d=inv(XDic'*XDic)*XDic'*y(1:N);
%% SPAMS  可持续加速算法
% d = mexConjGrad(XDic,y(1:N),dic,1e-3,100);
%%  共轭梯度下降
% a_ij=cell(1,1);
% a_ij{1,1}=sparse(x(:));
% epsilon=1e-6;
% [d] = conjgrad(a_ij,y(1:N),dic,epsilon);
end
function erro=Obj_min(y,XDic,dic)
    erro=norm(y-XDic*dic,'fro');
end



