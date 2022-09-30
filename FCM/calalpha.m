function [ alpha ] = calalpha( alpha,sdata,Best,g )
%更新衰减函数
[Nc,Nt,Ns] = size(sdata);%这里Nt+1，因为需要知道第0时刻的状态值
Nk = size(Best,3);%已经求得最优解的FCM个数
for j = 1:Nc
    if Nk == 0
        alpha(Nk+1,j) = 1;
    else
        Bestt = permute(Best(j,:,:),[3,2,1]);%横坐标是Nk,纵坐标是Nc
        for i = 1:Ns
            Alpha = repmat(alpha(1:Nk,j)',[Nt-1,1]);%Nt*Nk
            data = sdata(:,1:Nt-1,i);%Nc*Nt
            data = data';%Nt*Nc
            rdata = SFCM(data*Bestt',g);%Nt*Nk
            rrdata = sdata(j,2:Nt,i)' - sum(Alpha.*rdata,2);%Nt*1
            obj1 = sum(rrdata,1);
            alpha(Nk+1,j) = alpha(Nk+1,j) + obj1;
        end
        if alpha(Nk+1,j)>0
            alpha(Nk+1,j) = exp(-alpha(Nk+1,j));
        else
            alpha(Nk+1,j) = -exp(alpha(Nk+1,j));
        end
    end
end
end

