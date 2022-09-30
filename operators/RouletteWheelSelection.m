function index = RouletteWheelSelection(N,Fitness)
% RouletteWheelSelection���̶�ѡ����Ӧ��ֵԽ�󣬱�ѡ�����Խ��
% Input:Nѡ��N�����壬Fitness��Ӧ������
% Output:indexѡ�����±�
    Fitness = reshape(Fitness,1,[]);
    if (~isempty(find(Fitness<1, 1)))
        Fitness(Fitness == 0) = inf;
        Fitness = 1/min(Fitness)*Fitness;%��һ��
    end
    Fitness = cumsum(Fitness);
    Fitness = Fitness./max(Fitness);
    index   = arrayfun(@(S)find(rand<=Fitness,1),1:N);%�ҵ�rand<=Fitness�����ȳ��ֵ�1����Ϊ����±�
end