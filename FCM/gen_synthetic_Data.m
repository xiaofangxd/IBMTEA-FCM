function [ sFCM,sdata ] = gen_synthetic_Data( Ns,Nc,Nt,density,g )
% �˺������ںϳ�FCM
%   1,�������FCMȨ��
%   2.�����ʼ��ÿ���ڵ��Ӧ��״ֵ̬
%   3.����Ȩ�غͳ�ʼ��״ֵ̬������������
sFCM = sprand(Nc,Nc,density)+0;
sFCM(sFCM ~= 0) = sFCM( sFCM ~= 0 )*2 - 1;
while(sum(sum(abs(sFCM)<=0.05 & sFCM ~= 0)) ~= 0 )
    temp = rand(Nc,Nc)*2-1;
    sFCM(abs(sFCM)<=0.05 & sFCM ~= 0) = temp(abs(sFCM)<=0.05 & sFCM ~= 0);%��Ȩ��ֵ����ֵС��0.05����ʱ
end

sdata = rand(Nc,Nt+1,Ns);% ��һ�γ�ʼ�����ڶ��ΰ��յ�һ�β������Դ�����
for j = 1:Ns
    for i =2:Nt+1
        sdata(:,i,j) = SFCM(sFCM*sdata(:,i-1,j),g);
    end
end
end

