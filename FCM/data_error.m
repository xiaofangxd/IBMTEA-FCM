function [ obj] = data_error( x,adept,sdata,g,alpha,Best,batch)
% 此函数用来计算当前权重矩阵x对应的数据错误率
% x代表权重矩阵N*Nc
% sdata代表数据，行数代表FCM节点个数,列数代表时间序列长度
% g代表sigmoid函数的参数
% alpha是衰竭系数
% Best是前几个训练好的模糊认知图权重
[Nc,Nt,Ns] = size(sdata);%这里Nt+1，因为需要知道第0时刻的状态值
[N,Ncc] = size(x);
Nk = size(Best,3);%已经求得最优解的FCM个数
if Ncc ~= Nc
    disp("error!");
end
Bestt = permute(Best(adept,:,:),[3,2,1]);%横坐标是Nk,纵坐标是Nc
obj = zeros(N,1);
for i = 1:Ns
    Alpha = repmat(alpha(1:Nk,adept)',[Nt-1,1]);%Nt*Nk
    data = sdata(:,1:Nt-1,i);%Nc*Nt
    data = data';%Nt*Nc
    rdata = SFCM(data*Bestt',g);%Nt*Nk
    data = SFCM(data*x',g);%Nt*N
    if Nk ~= 0      
        rrdata = repmat(sdata(adept,2:Nt,i)' - sum(Alpha.*rdata,2),[1,N]);
%         rrrdata = abs(rrdata);
        rrrdata = rrdata;
        rrrdata(alpha(Nk+1,adept).*rrdata <= 0) = 0;
        data = alpha(Nk+1,adept).*data;
        obj2 = sum((data - rrrdata).^2,1);%最小化1*N
        obj = obj + obj2';
    else
        rrdata = repmat(sdata(adept,2:Nt,i)',[1,N]);
        obj2 = sum((data - rrdata).^2,1);%最小化1*N
        obj = obj + obj2';
    end
end
obj = obj./(Nc.*Ns.*(Nt-1));%N*1
end

