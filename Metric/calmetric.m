function [data_error,out_of_Sample_error,SS_Mean,model_error] = calmetric(Nk,ssdata,ALPHA,sFCM,var,Nc,Nt,Ns,g,flag)
% ����data_error,out_of_Sample_error,SS_Mean,model_error
% Input: ssdataѵ��ʱ���������ݡ�ALPHA˥��������sFCMs�Ǻϳɵ����ݡ�oBFCM��ѵ���õ�������ά��ΪNc*Nc*Nk��Nc�ǽڵ�����Nt�����ݳ��ȡ�Ns�ǳ�ʼ�������ĸ���
% Output: out_of_Sample_error,SS_Mean,model_error
%--------------------------------------------------------------------------
    data_error = 0;    
    out_of_Sample_error = 0;
    Nss = Ns;
    oBFCM = var;
    if flag == 1 %�ϳ�����
    %     Nk = size(ALPHA,1);
        sdata = rand(Nc,Nss);%��һ�γ�ʼ�����ڶ��ΰ��յ�һ�β������Դ�����
        sdata1 = sdata;
        for i=1:Nt
            sdata2 = permute(ssdata(:,i,:),[1,3,2]);
            odata = SFCM(sFCM*sdata,g);%Nc*Ns����sdata��ֵ����sFCM��Ӧ��ʱ������
            odata1 = zeros(Nc,Ns);
            odata2 = zeros(Nc,Ns);
            for j=1:Nk
                odata1 = odata1 + repmat(ALPHA(j,:)',1,Nss).*SFCM((oBFCM(:,:,j)*sdata1),g);%Nc*Nss����sdata��ֵ����oBFCM��Ӧ��ʱ������
                odata2 = odata2 + repmat(ALPHA(j,:)',1,Ns).*SFCM((oBFCM(:,:,j)*sdata2),g);%����ssdata��ֵ����oBFCM��Ӧ��ʱ������
            end
            sdata = odata;
            sdata1 = odata;
            data_error = data_error + sum(sum((permute(ssdata(:,i+1,:),[1,3,2]) - odata2).^2));
            out_of_Sample_error = out_of_Sample_error + sum(sum(abs(odata - odata1)));
        end
        data_error = data_error./(Ns*Nc*Nt);
        out_of_Sample_error = out_of_Sample_error./(Ns*Nc*Nt);
        if Nk == 1
            model_error = sum(sum(abs(sFCM-oBFCM)))./(Nc*Nc);
        else
            model_error = 0;
        end
        oBFCM(abs(oBFCM)>=0.05) = 1;%��Ȩ��ֵ����ֵ���ڵ���0.05ʱ����ֵΪ1
%         oBFCM = sum(oBFCM,3);
%         ooBFCM = oBFCM;
%         ooBFCM(oBFCM>=(Nk/2)) = 1;
%         ooBFCM(oBFCM<(Nk/2)) = 0;
        ooBFCM = oBFCM(:,:,1);
        sFCM(abs(sFCM)>=0.05) = 1;%��Ȩ��ֵ����ֵ���ڵ���0.05ʱ����ֵΪ1
        ntp = sum(sum(ooBFCM == sFCM & ooBFCM == 1));
        ntn = sum(sum(ooBFCM == sFCM & ooBFCM == 0));
        nfp = sum(sum(ooBFCM == 0 & sFCM == 1));
        nfn = sum(sum(ooBFCM == 1 & sFCM == 0));
        s1 = ntn./(ntn+nfp);s2 = ntp./(ntp+nfn);
        SS_Mean = 2*s1*s2/(s1+s2);
    elseif flag == 2 %Dream����
        model_error = 0;
        SS_Mean = 0;
        for i=1:Nt
            sdata2 = permute(ssdata(:,i,:),[1,3,2]);
            odata2 = zeros(Nc,Ns);
            for j=1:Nk
                odata2 = odata2 + repmat(ALPHA(j,:)',1,Ns).*SFCM((oBFCM(:,:,j)*sdata2),g);%����ssdata��ֵ����oBFCM��Ӧ��ʱ������
            end
            data_error = data_error + sum(sum((permute(ssdata(:,i+1,:),[1,3,2]) - odata2).^2));
        end
        data_error = data_error./(Ns*Nc*Nt);
    end
end