classdef TASK    
    %此类为所处理的任务信息，包含任务个数，任务维度，统一搜索空间，任务上下界和任务函数，此类需要用initTASK初始化
    properties
        M;%任务个数
        Tdims;%任务维度
        D_multitask;%统一搜索空间
        Lb;%任务的下界
        Ub;%任务的上界
        fun;%函数
    end    
    methods        
        function object = initTASK(object,Nc,sdata,g,alpha,Best,batch,bb)
            object.M = batch;%初始化任务个数
            object.Tdims = Nc*ones(object.M,1);%初始化任务维度
            object.D_multitask = Nc;%统一搜索空间
            object.Lb = -1*ones(object.M,object.D_multitask);%初始化下界
            object.Ub = ones(object.M,object.D_multitask);%初始化上界
            temp = (1+batch*(bb-1)):(bb*batch);
            for i = 1:object.M
                object.fun(i).fnc=@(x)data_error(x,temp(i),sdata,g,alpha,Best,batch);
            end
        end  
    end
end