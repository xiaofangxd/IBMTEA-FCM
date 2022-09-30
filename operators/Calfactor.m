function [objective,a,b] = Calfactor(objective,beta)
% ���������Ӧ�ȡ��������ء����صȼ�
% ���룺population��Ⱥ�ĸ�����Ϣ
% �����scalar_fitness������Ӧ�ȡ�skill_factor�������ء�factorial_ranks���صȼ���aÿ����������õ�ֵ��%��õ����ش��۶�Ӧ�ĸ���Ⱦɫ����Ϣ
    [N,~] = size(objective.factorial_costs);
    %�������صȼ�
    [x,y] = sort(objective.factorial_costs);%�������ش��۵�ֵ������������
    [~,objective.factorial_ranks] = sort(y);%������˳��Ϊ�����صȼ�
    %Population.factorial_ranks(y) = repmat([1:N]',[1,Task.M]);%������˳��Ϊ�����صȼ�
    a = x(1,:);%��õ����ش���
    bb = y(1,:)';%��õ����ش��۶�Ӧ�ĸ�������
    b = objective.rnvec(bb,:);%��õ����ش��۶�Ӧ�ĸ���Ⱦɫ����Ϣ
    
    %���㼼������
    [xx,yy] = min(objective.factorial_ranks,[],2);%yyΪ�����ó��������ţ�xxΪ�����ó�������ȼ�
    for i=1:N
        x = find(objective.factorial_ranks(i,:) == xx(i));%��i��Ⱦɫ��ı�����õ�����
        equivalent_skills = length(x);%��i��Ⱦɫ��ı�����õ�������
        if(equivalent_skills == 1)
            objective.skill_factor(i) = yy(i);
        else
            objective.skill_factor(i) = x(1+round((equivalent_skills-1)*rand(1)));
        end
        %���������Ӧ��
        objective.scalar_fitness(i) = 1./(beta.*objective.factorial_costs(i,objective.skill_factor(i))+1);
%         objective.scalar_fitness(i) = 1./(objective.factorial_costs(i,objective.skill_factor(i)));
    end
    %���������Ӧ��
        objective.scalar_fitness = 1./(beta*min(objective.factorial_ranks,[],2)+1);
%         objective.scalar_fitness = 1./(min(objective.factorial_ranks,[],2));
end