%%
function x0=LoCOMP(s,phi,k)
%% Code by Liu He ��aresmiki@163.com��
n=length(phi);
r0=s;
x0=zeros(size(s));
ln=length(r0);
%% ��һ���ڻ�
for j=1:ln-n+1
%         alp(j)=sum(r0(j:j+n-1).*phi);
    alp(j)=dot(r0(j:j+n-1),phi);
end

for i=1:k
     
    [~,indx(i)]=max(abs(alp)); %���ԭ��ѡ��
    val=alp(indx(i));
    if i>1
         sup=find(abs(indx(1:end-1)-indx(end))<n);  %֧��ԭ��
         if isempty(sup)
              x0(indx(i))=val;
              r0(indx(i):indx(i)+n-1)=r0(indx(i):indx(i)+n-1)-val*phi;
              %�����ڻ� 
              for j=max(indx(i)-n,1):min(indx(i)+n-1,ln-n+1)
                   alp(j)=dot(r0(j:j+n-1),phi);
              end
         else
           sup_des=sort([indx(sup),indx(end)]);
           L=(sup_des(end)+n-1-sup_des(1))+1;
           %%����ֲ�֧���ֵ�
           SupDic=get_Circshift(repmat([phi(:);zeros(L-n,1)],1,length(sup_des)),[sup_des-sup_des(1)]);
            x_temp=r0(sup_des(1):sup_des(end)+n-1);
            xcoff=pinv(SupDic)*x_temp(:);
            x0(sup_des)=x0(sup_des)+xcoff;
            r0(sup_des(1):sup_des(end)+n-1)=r0(sup_des(1):sup_des(end)+n-1)-SupDic*xcoff;
             %�����ڻ� 
            for j=max(sup_des(1)-n,1):min(sup_des(end)+n-1,ln-n+1)
                alp(j)=dot(r0(j:j+n-1),phi);
            end
         end
    else
        x0(indx(i))=val;
        r0(indx(i):indx(i)+n-1)=r0(indx(i):indx(i)+n-1)-val*phi;
        %�����ڻ� 
        for j=max(indx(i)-n,1):min(indx(i)+n-1,ln-n+1)
            alp(j)=dot(r0(j:j+n-1),phi);
        end
    end

end