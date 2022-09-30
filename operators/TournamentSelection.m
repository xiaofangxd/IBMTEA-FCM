function index = TournamentSelection(K,N,Fitness)
% TournamentSelection���ƽ�����(��Ӧ��ֵԽ�󣬱�ѡ����Խ��)
% Input: K��Ԫ��������Nѡ��N�����壬Fitness��Ӧ������
% Output��indexѡ�����������
    Fitness = reshape(Fitness,1,[]);
    if (~isempty(find(Fitness<1, 1)))
        Fitness(Fitness == 0) = inf;
        Fitness = 1/min(Fitness)*Fitness;%��һ��
    end
    [~,rank] = sort(Fitness,'descend');
    [~,rank] = sort(rank);
    Parents  = randi(length(Fitness),K,N);
    [~,best] = min(rank(Parents),[],1);
    index    = Parents(best+(0:N-1)*K);
end