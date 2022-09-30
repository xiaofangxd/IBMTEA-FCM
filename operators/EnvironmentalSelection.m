function nPopulation = EnvironmentalSelection(Population,selection_process,N,M)
% ����ѡ��
% Input:Population��Ⱥ��selection_processѡ�����(��ȡelitist��roulette
% wheel),Nѡ���������Ŀ,M������
% Output:nPopulationѡ��֮�����Ⱥ
%--------------------------------------------------------------------------
    nPopulation = INDIVIDUAL();
    if strcmp(selection_process,'elitist')%��Ӣѡ��
        [~,y]=sort(Population.scalar_fitness,'descend');
    elseif strcmp(selection_process,'roulette wheel')%���̶�ѡ��
        y = [];
        for i = 1:M
            ind = find(Population.skill_factor == i);%���ó��ļ��ܽ��з���ѡ��
            Fitness = Population.scalar_fitness(Population.skill_factor == i,:);%�ҵ������Ӧ�ĸ������Ӧ�Ⱥ���ֵ
            ind1 = ind(RouletteWheelSelection(N/M,Fitness));%ÿ������ѡ��N/M���壬�˴�ѡ���ĸ�������
            y = [y;ind1];
        end
    elseif strcmp(selection_process,'Tournament')%������ѡ��
        y = [];
        for i = 1:M
            ind = find(Population.skill_factor == i);%���ó��ļ��ܽ��з���ѡ��
            Fitness = Population.scalar_fitness(Population.skill_factor == i,:);%�ҵ������Ӧ�ĸ������Ӧ�Ⱥ���ֵ
            ind1 = ind(TournamentSelection(2,N/M,Fitness));%ÿ������ѡ��N/M���壬�˴�ѡ���ĸ���������Ĭ�϶�Ԫ������
            y = [y;ind1];
        end
    end
    nPopulation.rnvec = Population.rnvec(y(1:N),:);
    nPopulation.factorial_costs = Population.factorial_costs(y(1:N),:);
    nPopulation.factorial_ranks = Population.factorial_ranks(y(1:N),:); 
    nPopulation.scalar_fitness = Population.scalar_fitness(y(1:N),:);
    nPopulation.skill_factor = Population.skill_factor(y(1:N),:);
end