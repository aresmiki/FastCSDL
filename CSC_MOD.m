function [dic,xcofe,rec]=CSC_MOD(y,np,kp,iter,disp)
%% Code by Liu He ��aresmiki@163.com��

if (nargin <= 3)   %
    iter=100;
end
if (nargin <= 4)   %
    disp=2;
end
%% ���پ��ϡ�����
dic=randn(np,1);
dic=dic./norm(dic);
N=length(y(:));
erro(1)=0;

for i=1:iter     
        %% �̶�ԭ�ӣ�ϡ�����
        %% LoCOMP
        xcofe=LoCOMP(y,dic,kp);
        rec=Rec_LoCOMP(dic,xcofe,N);
    %% �����㣬��ֹ����    
        erro(i+1)=norm(y(:)-rec(:),'fro');
        if disp==1
          fprintf('iter = %d\n   erro = %f\n\r',i,erro(i+1));
        end
        if abs(erro(i+1)-erro(i))<1e-4
            break;
        end
        %% �̶�ϡ�裬����ԭ��
        dic=update_dic(y,xcofe,dic);
        dic=dic./norm(dic); 
end
