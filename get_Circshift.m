
function CRecq=get_Circshift(Recq,lambda)
%% Code by Liu He £¨aresmiki@163.com£©
CRecq=(Recq);
for i=1:length(lambda)
    CRecq(:,i) = circshift(Recq(:,i),lambda(i));
end




