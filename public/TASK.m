classdef TASK    
    %����Ϊ�������������Ϣ�������������������ά�ȣ�ͳһ�����ռ䣬�������½����������������Ҫ��initTASK��ʼ��
    properties
        M;%�������
        Tdims;%����ά��
        D_multitask;%ͳһ�����ռ�
        Lb;%������½�
        Ub;%������Ͻ�
        fun;%����
    end    
    methods        
        function object = initTASK(object,Nc,sdata,g,alpha,Best,batch,bb)
            object.M = batch;%��ʼ���������
            object.Tdims = Nc*ones(object.M,1);%��ʼ������ά��
            object.D_multitask = Nc;%ͳһ�����ռ�
            object.Lb = -1*ones(object.M,object.D_multitask);%��ʼ���½�
            object.Ub = ones(object.M,object.D_multitask);%��ʼ���Ͻ�
            temp = (1+batch*(bb-1)):(bb*batch);
            for i = 1:object.M
                object.fun(i).fnc=@(x)data_error(x,temp(i),sdata,g,alpha,Best,batch);
            end
        end  
    end
end