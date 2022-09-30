function [ obj] = data_error( x,adept,sdata,g,alpha,Best,batch)
% �˺����������㵱ǰȨ�ؾ���x��Ӧ�����ݴ�����
% x����Ȩ�ؾ���N*Nc
% sdata�������ݣ���������FCM�ڵ����,��������ʱ�����г���
% g����sigmoid�����Ĳ���
% alpha��˥��ϵ��
% Best��ǰ����ѵ���õ�ģ����֪ͼȨ��
[Nc,Nt,Ns] = size(sdata);%����Nt+1����Ϊ��Ҫ֪����0ʱ�̵�״ֵ̬
[N,Ncc] = size(x);
Nk = size(Best,3);%�Ѿ�������Ž��FCM����
if Ncc ~= Nc
    disp("error!");
end
Bestt = permute(Best(adept,:,:),[3,2,1]);%��������Nk,��������Nc
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
        obj2 = sum((data - rrrdata).^2,1);%��С��1*N
        obj = obj + obj2';
    else
        rrdata = repmat(sdata(adept,2:Nt,i)',[1,N]);
        obj2 = sum((data - rrdata).^2,1);%��С��1*N
        obj = obj + obj2';
    end
end
obj = obj./(Nc.*Ns.*(Nt-1));%N*1
end

