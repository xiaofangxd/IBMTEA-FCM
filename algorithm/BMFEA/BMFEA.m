function data_BMFEA = BMFEA(rmp,Pc,Pm,Pd,mu,muu,sigma,selection_process,beta,probswap)
    % 本程序主要实现了用进化多任务优化算法优化多个模糊认知图参数、优化函数为最小化函数、最大化函数需要转化为最小化函数
    % Input:  Nk是FCM个数，rmp随机交配池概率,pi_l个体学习的概率,Pc模拟二进制交叉概率,mu模拟二进制交叉参数(可调),Pm多项式变异概率,muu多项式变异参数(可调),sigma高斯变异模型的标准差(可调),
    %         selection_process可供选择：elitist、roulette wheel、Tournament,
    %         name构建FCM时可选择的转换函数有：SFCM、WFCM、HTFCM、TLFCM
    %         options调用matlab函数优化器（拟牛顿法）->设置预学习优化器
    %         beta：data_error参数、probswap：位交换概率
    % Output: data_MFEA（运行时间、每代最好值、最好个体索引（调试使用的，因为每代都变，所以这个参数没啥意义）、总评价次数）
    global Nk N gen Nc sdata g batch
    tic
    Ncc = Nc/batch;
    NN=0;
    %% 0.衰减系数矩阵、最优个体矩阵
    alpha = zeros(Nk,Nc);
    Best = zeros(Nc,Nc,Nk);
    EvBestFitness = zeros(gen+1,batch,Nk);           %每代最好的适应度值
    TotalEvaluations=zeros(gen+1,1,Nk);              %每代每个个体评价次数
    allEvBestFitness = zeros(gen+1,Nc,Nk);           %所有任务每代最好的适应度值
    %% 1.分批次学习Nk个FCM
    for j = 1:Nk
        % 1.1 初始化FCM和任务并计算alpha
        disp(['Nk = ',num2str(j)]);
        alpha = calalpha(alpha,sdata,Best(:,:,1:j-1),g);
        for bb = 1:Ncc       
            disp(['batch = ',num2str(bb)]);
            Task = TASK();
            Task = initTASK(Task,Nc,sdata,g,alpha(1:j,:),Best(:,:,1:j-1),batch,bb);

            % 1.2 初始化种群
            Population = INDIVIDUAL();                    %生成初始种群
            Population = initPOP(Population,N,Task.D_multitask,Task.M);

            % 1.3 根据多任务环境中的每个优化任务评估每个个体的因子代价
            [Population,TotalEvaluations(1,1,j)] = evaluate(Population,Task);

            % 1.4 计算初始化种群的因素等级以及技能因素
            [Population,EvBestFitness(1,:,j),bestind] = Calfactor(Population,beta);

            % 1.5 优化过程
            for i = 1:gen
                %4.1 个体变异交叉
                Offspring  = GA_MFEA(Population,rmp,Pc,Pm,Pd,mu,muu,sigma,probswap,i,gen);
                %4.2 计算因子代价
                [Offspring,TotalEvaluations(i+1,1,j)] = evaluate(Offspring,Task);
                TotalEvaluations(i+1,1,j) = TotalEvaluations(i+1,1,j) + TotalEvaluations(i+1,1,j);
                %4.3 种群合并
                intpopulation = combpop(Population,Offspring);
                %4.4 更新标量适应度，技能因素，因素等级
                [intpopulation,EvBestFitness(i+1,:,j),bestind] = Calfactor(intpopulation,beta);
                %4.5 环境选择
                Population = EnvironmentalSelection(intpopulation,selection_process,N,Task.M);
                disp(['BMFEA Generation = ', num2str(i), ' EvBestFitness = ', num2str(EvBestFitness(i+1,:,j))]);%为了记录初始化的值所以次数+1
            end

            % 1.6更新模糊认知图的评价方式
            Best((1+batch*(bb-1)):(bb*batch),:,j) = bestind;
            allEvBestFitness(:,(1+batch*(bb-1)):(bb*batch),:) = EvBestFitness;
        end
        % 1.7满足条件退出
        if j ~= 1
            temp = alpha;temp(alpha<0) = log(-temp(alpha<0));temp(alpha>0) = -log(temp(alpha>0));
            temp1 = sum(temp(j,:));temp2 = sum(temp(j-1,:));
            if abs(temp1-temp2)<=1e-3
                NN = j;
                break;
            else
                NN = Nk;
            end
        else
            NN = Nk;
        end
    end

    %% 2.记录算法结果
    data_BMFEA.Nk = NN;
    data_BMFEA.wall_clock_time=toc;
    data_BMFEA.alpha=alpha;
    data_BMFEA.EvBestFitness=EvBestFitness;
    data_BMFEA.allEvBestFitness=allEvBestFitness;
    data_BMFEA.bestInd_data=Best;
    data_BMFEA.TotalEvaluations=TotalEvaluations;
end