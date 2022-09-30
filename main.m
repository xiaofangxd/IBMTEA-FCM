% ��������Ҫʵ���˽����������Ż��㷨�Ż�����ģ����֪ͼ��Ȩ�ء��Ż�����Ϊ��С�������ʺ���
% ���κ����������ϵ�ҵ�����: Xiao Feng(Email: xiaofengxd@126.com��
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
    %% FCM��������
    Nc = 500;% FCM�Ľڵ����,20,30
    Nk = 1;%FCM�ĸ�������ѧϰ
    density = 0.4; % �ϳɵ�FCM�з���Ȩ����ռ����
    Ns = 1;% ��ʼ�������ĸ���
    Nt = 20;% �ϳ����ݵĸ���
    g = 0.5;% SFCM��ת�������Ĳ�����Nc>=200,g=0.5��
    beta = 1000;% ���������ʱ�Ĳ�����Nc>=200,beta = 1000
    flag = 1;% 1�ܺϳ����ݣ�2��DREAM

    %% BIMTEA��������
    gen = 200;                                  % ��������
    batchsize = [5,1,Nc];
    batch = batchsize(ppp);                     % ÿ��ѧϰFCM�ĸ���
    N = 100;                                    % ��Ⱥ��С������Ϊż����
    rmp = 0.3;                                  % �������ظ���
    Pc = 1;                                     % ģ������ƽ������
    Pm = 1;                                     % ����ʽ�������
    Pd = 0.5;                                   % dropout����
    mu = 20;                                    % ģ������ƽ������(�ɵ�)
    muu = 20;                                   % ����ʽ�������
    probswap = 0.5;                             %������������
    sigma = 0.02;                               % ��˹����ģ�͵ı�׼��(�ɵ�)   
    selection_process = 'elitist';              % �ɹ�ѡ��elitist��roulette wheel��Tournament
    times = 5;                                  %�㷨���д���
    
    allEvBestFitness = zeros(gen+1,Nc,Nk);           %��ͼ��������ÿ����õ���Ӧ��ֵ

    %% �����ϳ�����
    [sFCM, sdata] = gen_synthetic_Data(Ns,Nc,Nt,density,g);

    %% ����DREAM3��DREAM4
    if flag == 2
        Nt = 10;%���������
        Dreamnum = 4;%���ݼ���ţ���ѡ3��4
        Nc = 100;%�ڵ�����Dream3��ѡ10,50,100,Dream4��ѡ10,100
        datanum = 5;%ÿ�����ݼ������һ�飬��1��5
        Ns = 10;%DREAM3��Nc=10��ӦNs=4,DREAM3��Nc=50��ӦNs=23,DREAM3��Nc=100��ӦNs=46,DREAM4��Nc=10��ӦNs=5,DREAM4��Nc=100��ӦNs=10
        fid = fopen(['Dream\D',num2str(Dreamnum),'_',num2str(Nc),'_',num2str(datanum),'.tsv']);
        FormatString=repmat('%f ',1,Nc+1);
        out =cell2mat(textscan(fid,FormatString,Ns*21,'HeaderLines',1));
        sdata = reshape(out(:,2:end)',[Nc,21,size(out(:,2:end)',2)/21]);
        sdata = sdata(:,11:end,:);
    end

    %% ���ʵ��
    out_of_Sample_error = zeros(times,1);SS_Mean = zeros(times,1);data_error = zeros(times,1);model_error = zeros(times,1);
    for i =1:times
        disp(['times = ',num2str(i)]);
        data_BMFEA(i) = BMFEA(rmp,Pc,Pm,Pd,mu,muu,sigma,selection_process,beta,probswap);%����BMFEA�㷨
        allEvBestFitness = allEvBestFitness + data_BMFEA(i).allEvBestFitness;
        vars = data_BMFEA(i).bestInd_data;
        ALPHA = data_BMFEA(i).alpha;
        disp(['Nk=',num2str(data_BMFEA(i).Nk)]);
        [data_error(i),out_of_Sample_error(i),SS_Mean(i),model_error(i)] = calmetric(data_BMFEA(i).Nk,sdata,ALPHA,sFCM,vars,Nc,Nt,Ns,g,flag);
    end