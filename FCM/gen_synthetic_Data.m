function [ sFCM,sdata ] = gen_synthetic_Data( Ns,Nc,Nt,density,g )
% 此函数用于合成FCM
%   1,随机构造FCM权重
%   2.随机初始化每个节点对应的状态值
%   3.根据权重和初始化状态值产生数据序列
sFCM = sprand(Nc,Nc,density)+0;
sFCM(sFCM ~= 0) = sFCM( sFCM ~= 0 )*2 - 1;
while(sum(sum(abs(sFCM)<=0.05 & sFCM ~= 0)) ~= 0 )
    temp = rand(Nc,Nc)*2-1;
    sFCM(abs(sFCM)<=0.05 & sFCM ~= 0) = temp(abs(sFCM)<=0.05 & sFCM ~= 0);%当权重值绝对值小于0.05等于时
end

sdata = rand(Nc,Nt+1,Ns);% 第一次初始化，第二次按照第一次产生，以此类推
for j = 1:Ns
    for i =2:Nt+1
        sdata(:,i,j) = SFCM(sFCM*sdata(:,i-1,j),g);
    end
end
end

