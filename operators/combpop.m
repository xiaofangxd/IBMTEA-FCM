function [ intpopulation ] = combpop( Population,Offspring )
%�ϲ��������Ӵ���Ⱥ
% input:������Ⱥ���Ӵ���Ⱥ
% output:��������Ⱥ
intpopulation = INDIVIDUAL();
intpopulation.rnvec = [Population.rnvec;Offspring.rnvec];
intpopulation.factorial_costs = [Population.factorial_costs;Offspring.factorial_costs];
intpopulation.factorial_ranks = [Population.factorial_ranks;Offspring.factorial_ranks];
intpopulation.skill_factor = [Population.skill_factor;Offspring.skill_factor];
intpopulation.scalar_fitness = [Population.scalar_fitness;Offspring.scalar_fitness];
end

