function Offspring = GA_MFEA(Parent,rmp,Pc,Pm,Pd,disC,disM,sigma,probswap,k,gen)
% 此函数功能是通过模拟二进制交叉和高斯变异产生子代，并利用垂直文化传播进行技能因子的继承（两两组成的父代个体，必须进行交叉或者变异）。
% Input: Parent父代信息（染色体，技能因子）、rmp文化交流参数、Pc模拟二进制交叉概率、disC交叉参数、sigma高斯变异参数
% Output: Offspring子代信息（染色体，技能因子）
% 第“1.”模式中，通过变异产生的两个后代可能具有相同的父母；第“2.”模式中，通过变异产生的两个后代父母一定不同
    [N,~] = size(Parent.rnvec);
    select = randperm(N);
    rrnvec = Parent.rnvec(select,:);%打乱顺序
    sskill_factor = Parent.skill_factor(select,:); 
    Parent1 = rrnvec(1:floor(end/2),:);
    factor1 = sskill_factor(1:floor(end/2),:);
    Parent2 = rrnvec(floor(end/2)+1:floor(end/2)*2,:);
    factor2 = sskill_factor(floor(end/2)+1:floor(end/2)*2,:);
    Offspring = INDIVIDUAL();
    %Offspring.skill_factor = zeros(N,1);%1.初始化子代的技能因子为0
    Offspring.skill_factor = sskill_factor;%2.初始化子代的技能因子对应为父代的技能因子
    factorb1 = repmat(1:N/2,1,2);
    factorb2 = repmat(N/2+1:N,1,2);
    temp = randi(2,1,N);%对于子代随机选择它是继承第一个父母还是第二个父母
    offactor = zeros(1,N);
    offactor(temp == 1) = factorb1(temp == 1);
    offactor(temp == 2) = factorb2(temp == 2);%子代继承父母的编号
    %Offspring.skill_factor = sskill_factor(offactor);%1.所有子代继承父代基因
    [NN,D]   = size(Parent1);
    
     %% 对于factor1 == factor2 or rand<RMP
    % Simulated binary crossover
    beta = zeros(NN,D);
    mu   = rand(NN,D);
    rmpp = rand(NN,1);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta(repmat(factor1 ~= factor2 & rmpp>=rmp,1,D)) = 1;%不同技能因子的个体满足rmp交叉,在模拟二进制交叉的时候beta=1为不交叉
    beta(repmat(rand(NN,1)>=Pc,1,D)) = 1;
    Offspring.rnvec = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                 (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
    Offspring.skill_factor(repmat(beta(:,1) ~= 1,2,1)) = sskill_factor(offactor(repmat(beta(:,1) ~= 1,2,1)));%2.对于交叉的个体，其随机选择一个父代继承基因。
    
    % Polynomial mutation
    Lower = -ones(2*NN,D);
    Upper = ones(2*NN,D);
    Site  = rand(2*NN,D) < Pm/D;
    mu    = rand(2*NN,D);
    temp  = Site & mu<=0.5 & repmat(beta(:,1) ~= 1,2,D);
    Offspring.rnvec = min(max(Offspring.rnvec,Lower),Upper);
    Offspring.rnvec(temp) = Offspring.rnvec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspring.rnvec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5 & repmat(beta(:,1) ~= 1,2,D); 
    Offspring.rnvec(temp) = Offspring.rnvec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspring.rnvec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
                
%     % nonuniform mutation
%     temp2 = rand(NN*2,D);
%     yita = 1 - temp2.^((1-k/gen).^5);
%     temp3 = rand(NN*2,D);temp3(temp3>=0.5)=1;temp3(temp3<0.5)=0;
%     Offspring.rnvec(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 0) = Offspring.rnvec(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 0) + (1 - Offspring.rnvec(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 0)).*yita(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 0);%2.只对没有交叉的个体进行变异
%     Offspring.rnvec(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 1) = Offspring.rnvec(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 1) - (1 + Offspring.rnvec(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 1)).*yita(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 1);

    
    % variable swap (uniform X)
    swap_indicator = (rand(NN,D) >= probswap);
    temp = swap_indicator & repmat(factor1 == factor2,1,D);
    temp1 = Offspring.rnvec(1:NN,:);temp2 = Offspring.rnvec(NN+1:end,:);
    temp3 = temp2;
    temp3(temp) = temp1(temp);
    temp1(temp) = temp2(temp);
    Offspring.rnvec(1:NN,:) = temp1;Offspring.rnvec(NN+1:end,:) = temp3;
    
    %% 对于factor1 ~= factor2 and rand>=RMP
    rrrnvec = rrnvec;
%     tempp = beta(:,1) == 1;
    tempp = factor1 ~= factor2 & rmpp>=rmp;
    aa = sskill_factor(tempp);
    bb = find(tempp);
    for i = 1:length(aa)
        cc = find(sskill_factor == aa(i));
        if length(cc) ~= 1
            cc(cc == bb(i)) = [];
        end
        rrrnvec(bb(i),:) = rrnvec(randi(length(cc),1),:);
    end
    % Simulated binary crossover
    beta = zeros(N,D);
    mu   = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta(repmat(~tempp,2,D)) = 1;
    beta(repmat(rand(N,1)>=Pc,1,D)) = 1;
    Offspringrnvec = [(rrnvec+rrrnvec)/2+beta.*(rrnvec-rrrnvec)/2
                 (rrnvec+rrrnvec)/2-beta.*(rrnvec-rrrnvec)/2];
    
    % Polynomial mutation
    Lower = -ones(2*N,D);
    Upper = ones(2*N,D);
    Site  = rand(2*N,D) < Pm/D;
    mu    = rand(2*N,D);
    temp  = Site & mu<=0.5 & repmat(beta(:,1) ~= 1,2,D);
    Offspringrnvec = min(max(Offspringrnvec,Lower),Upper);
    Offspringrnvec(temp) = Offspringrnvec(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                      (1-(Offspringrnvec(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5 & repmat(beta(:,1) ~= 1,2,D); 
    Offspringrnvec(temp) = Offspringrnvec(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                      (1-(Upper(temp)-Offspringrnvec(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
                  
%     % nonuniform mutation
%     temp2 = rand(N*2,D);
%     yita = 1 - temp2.^((1-k/gen).^5);
%     temp3 = rand(N*2,D);temp3(temp3>=0.5)=1;temp3(temp3<0.5)=0;
%     Offspringrnvec(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 0) = Offspringrnvec(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 0) + (1 - Offspringrnvec(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 0)).*yita(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 0);
%     Offspringrnvec(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 1) = Offspringrnvec(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 1) - (1 + Offspringrnvec(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 1)).*yita(repmat(beta(:,1) ~= 1,2,D) & temp2 <= Pm & temp3 == 1);
    
    % variable swap (uniform X)
    swap_indicator = (rand(N,D) >= probswap);
    temp = swap_indicator & repmat(beta(:,1) ~= 1,1,D);
    temp1 = Offspringrnvec(1:N,:);temp2 = Offspringrnvec(N+1:end,:);
    temp3 = temp2;
    temp3(temp) = temp1(temp);
    temp1(temp) = temp2(temp);
    Offspringrnvec(1:N,:) = temp1;Offspringrnvec(N+1:end,:) = temp3;
    Offspring.rnvec(bb,:) = Offspringrnvec(bb,:);
    
%     % Norm mutation
%     rvec=normrnd(0,sigma,[N,D]);
%     Offspring.rnvec(repmat(beta(:,1) == 1,2,D)) = Offspring.rnvec(repmat(beta(:,1) == 1,2,D)) + rvec(repmat(beta(:,1) == 1,2,D));%2.只对没有交叉的个体进行变异    
    Offspring.rnvec(Offspring.rnvec>1)=1;
    Offspring.rnvec(Offspring.rnvec<-1)=-1;
    Offspring.rnvec(abs(Offspring.rnvec)<0.05)=0;
    
    % dropout正则化
    Offspring.rnvec(rand(N,D)<Pd) = 0;
end