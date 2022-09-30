% 本程序主要实现了进化多任务优化算法优化集成模糊认知图的权重、优化函数为最小化错误率函数
% 有任何问题可以联系我的邮箱: Xiao Feng(Email: xiaofengxd@126.com）
clc,clear all
global Nk N gen Nc sdata g batch
% Defaults for this blog post
width = 3;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
% for ppp = 1:3
    ppp = 3;
    %% FCM参数设置
    Nc = 500;% FCM的节点个数,20,30
    Nk = 1;%FCM的个数集成学习
    density = 0.4; % 合成的FCM中非零权重所占比例
    Ns = 1;% 初始化向量的个数
    Nt = 20;% 合成数据的个数
    g = 0.5;% SFCM中转换函数的参数，Nc>=200,g=0.5；
    beta = 1000;% 计算错误率时的参数，Nc>=200,beta = 1000
    flag = 1;% 1跑合成数据，2跑DREAM

    %% BIMTEA参数设置
    gen = 200;                                  % 迭代次数
    batchsize = [5,1,Nc];
    batch = batchsize(ppp);                     % 每批学习FCM的个数
    N = 100;                                    % 种群大小（设置为偶数）
    rmp = 0.3;                                  % 随机交配池概率
    Pc = 1;                                     % 模拟二进制交叉概率
    Pm = 1;                                     % 多项式变异概率
    Pd = 0.5;                                   % dropout正则化
    mu = 20;                                    % 模拟二进制交叉参数(可调)
    muu = 20;                                   % 多项式交叉参数
    probswap = 0.5;                             %变量交换概率
    sigma = 0.02;                               % 高斯变异模型的标准差(可调)   
    selection_process = 'elitist';              % 可供选择：elitist、roulette wheel、Tournament
    times = 5;                                  %算法运行次数
    
    allEvBestFitness = zeros(gen+1,Nc,Nk);           %画图所有任务每代最好的适应度值

    %% 产生合成数据
    [sFCM, sdata] = gen_synthetic_Data(Ns,Nc,Nt,density,g);

    %% 读入DREAM3和DREAM4
    if flag == 2
        Nt = 10;%不允许更改
        Dreamnum = 4;%数据集编号，可选3，4
        Nc = 100;%节点数，Dream3可选10,50,100,Dream4可选10,100
        datanum = 5;%每个数据集里的哪一组，从1到5
        Ns = 10;%DREAM3的Nc=10对应Ns=4,DREAM3的Nc=50对应Ns=23,DREAM3的Nc=100对应Ns=46,DREAM4的Nc=10对应Ns=5,DREAM4的Nc=100对应Ns=10
        fid = fopen(['Dream\D',num2str(Dreamnum),'_',num2str(Nc),'_',num2str(datanum),'.tsv']);
        FormatString=repmat('%f ',1,Nc+1);
        out =cell2mat(textscan(fid,FormatString,Ns*21,'HeaderLines',1));
        sdata = reshape(out(:,2:end)',[Nc,21,size(out(:,2:end)',2)/21]);
        sdata = sdata(:,11:end,:);
    end

    %% 多次实验
    out_of_Sample_error = zeros(times,1);SS_Mean = zeros(times,1);data_error = zeros(times,1);model_error = zeros(times,1);
    for i =1:times
        disp(['times = ',num2str(i)]);
        data_BMFEA(i) = BMFEA(rmp,Pc,Pm,Pd,mu,muu,sigma,selection_process,beta,probswap);%调用BMFEA算法
        allEvBestFitness = allEvBestFitness + data_BMFEA(i).allEvBestFitness;
        vars = data_BMFEA(i).bestInd_data;
        ALPHA = data_BMFEA(i).alpha;
        disp(['Nk=',num2str(data_BMFEA(i).Nk)]);
        [data_error(i),out_of_Sample_error(i),SS_Mean(i),model_error(i)] = calmetric(data_BMFEA(i).Nk,sdata,ALPHA,sFCM,vars,Nc,Nt,Ns,g,flag);
    end