function d=update_dic(y,x,dic)
%% Code by Liu He ��aresmiki@163.com��
% x  ��ϵ���� y ����Ҫ���Ƶĺ��� �� np ���ֵ�ԭ�Ӵ�С
%% �� x�������ֵ�
np=length(dic);
N=size(x(:),1);
%%
% XDic=circulant(x(:));
% XDic=XDic(:,1:np);
XDic=get_Circshift(repmat(x(:),1,np),[0:1:np-1]); 
y=y(:);
%  XDic*dic �Ǿ�����. %�ⲿ�ֿ�ֱ���þ���������˷�
%% argmin ||y-XDic*dic||_F 
% %% ���ݶ��½����
% options = optimoptions(@fmincon,'HessianApproximation','lbfgs','MaxIterations',20);
% d=fmincon(@(xv)Obj_min(y(1:N),XDic,xv),dic,[],[],[],[],[],[],[],options);
% d=cgls(XDic,y(1:N),0,1e-3,20,'true',dic);
d=cgls(XDic,y(1:N),0,1e-3,20,[],dic); %�����ݶ��½�
% d=inv(XDic'*XDic)*XDic'*y(1:N);
%% SPAMS  �ɳ��������㷨
% d = mexConjGrad(XDic,y(1:N),dic,1e-3,100);
%%  �����ݶ��½�
% a_ij=cell(1,1);
% a_ij{1,1}=sparse(x(:));
% epsilon=1e-6;
% [d] = conjgrad(a_ij,y(1:N),dic,epsilon);
end
function erro=Obj_min(y,XDic,dic)
    erro=norm(y-XDic*dic,'fro');
end



