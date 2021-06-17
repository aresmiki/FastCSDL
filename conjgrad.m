% Convolutional Conjugate Gradient Least Squares
% E. Plaut, R. Giryes, "A Greedy Approach to Convolutional Sparse Coding", 2018

function [d] = conjgrad(a_ij,p,D,epsilon)
    d=D(:);
    [m2,n]=size(D);
    m=sqrt(m2);
    M=size(a_ij{1,1},1)-m+1;
    N=size(a_ij{1,1},2)-m+1;
    s_i=size(a_ij,1);
    Qd=zeros(m2,n);
    for i=1:s_i
        Aconvs=cell(n,1);
        a_flip=cell(n,1);
        for j=1:n
            Aconvs{j,1}=conv2(reshape(D(:,j),[m m]),reshape(full(a_ij{i,j}),[M+m-1,N+m-1]));
            a_flip{j,1}=flip(flip(reshape(a_ij{i,j},[M+m-1,N+m-1]),2));
        end
        for j1=1:n
            for j2=1:n
                flipAconvAconvs=conv2(Aconvs{j2,1},full(a_flip{j1,1}));
                flipAconvAconvs=flipAconvAconvs(M+m-1:M+m-1+m-1,N+m-1:N+m-1+m-1);
                Qd(:,j1)=Qd(:,j1)+flipAconvAconvs(:);
            end
        end
    end
    Qd=Qd(:);
    r = p - Qd;
    s = r;
    rsold = r' * r;
    for iter = 1:length(p)
        Qs=zeros(m2,n);
        for i=1:s_i
            Aconvs=cell(n,1);
            a_flip=cell(n,1);
            for j=1:n
               Aconvs{j,1}=conv2(reshape(s(1+m2*(j-1):m2*j),[m m]),reshape(full(a_ij{i,j}),[M+m-1,N+m-1]));
               a_flip{j,1}=flip(flip(reshape(a_ij{i,j},[M+m-1,N+m-1]),2));
            end
            for j1=1:n
                for j2=1:n
                    flipAconvAconvs=conv2(Aconvs{j2,1},full(a_flip{j1,1}));
                    flipAconvAconvs=flipAconvAconvs(M+m-1:M+m-1+m-1,N+m-1:N+m-1+m-1);
                    Qs(:,j1)=Qs(:,j1)+flipAconvAconvs(:);
                end
            end
        end
        Qs=Qs(:);
        alpha = rsold / max((s' * Qs),1e-4);
        d = d + alpha * s;
        r = r - alpha * Qs;
        rsnew = r' * r;
        if rsnew < epsilon
              break;
        end
        s = r + (rsnew / rsold) * s;
        rsold = rsnew;
    end
end