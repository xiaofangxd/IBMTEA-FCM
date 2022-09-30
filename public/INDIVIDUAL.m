classdef INDIVIDUAL  
    %�������һ����Ⱥ�������������ɣ�Ⱦɫ�塢���ش��ۡ����صȼ���������Ӧ�ȡ��������أ��ᡪ��Ⱥ��С���ݡ�ά��/������
    %��Ⱥ��Ҫ��initPOP��ʼ����������Ⱥ����ʱ��Ҫevaluate����
    properties
        rnvec; % (genotype)--> decode to find design variables --> (phenotype) 
        factorial_costs;%���ش���
        factorial_ranks;%���صȼ�
        scalar_fitness;%������Ӧ��
        skill_factor;%��������
    end    
    methods        
        function object = initPOP(object,N,D,MM)  
%             object.rnvec = rand(N,D);%��ʼ���������
            object.rnvec = rand(N,D)*2-1;%��ʼ���������
            object.rnvec(abs(object.rnvec)<0.05) = 0;%��Ȩ��ֵ����ֵС��0.05ʱ����ֵΪ0
            object.factorial_costs = inf*ones(N,MM);%��ʼ���������ش���Ϊ0
            object.factorial_ranks = zeros(N,MM);%��ʼ���������صȼ�Ϊ0
            object.scalar_fitness = zeros(N,1);%��ʼ�����������Ӧ��Ϊ0
            object.skill_factor = zeros(N,1);%��ʼ������ļ�������Ϊ0
        end
        
        function [object,call] = evaluate(object,Task)%��Ӧ������
            object.factorial_costs(:)=inf;
            call = 0;
            for i = 1:Task.M
                [object.factorial_costs(:,i),object.rnvec,calls]=CalObj(Task,object.rnvec,i,object.skill_factor);
                call = call + calls;
            end
        end
    end
end