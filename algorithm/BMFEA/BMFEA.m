function data_BMFEA = BMFEA(rmp,Pc,Pm,Pd,mu,muu,sigma,selection_process,beta,probswap)
    % ��������Ҫʵ�����ý����������Ż��㷨�Ż����ģ����֪ͼ�������Ż�����Ϊ��С����������󻯺�����Ҫת��Ϊ��С������
    % Input:  Nk��FCM������rmp�������ظ���,pi_l����ѧϰ�ĸ���,Pcģ������ƽ������,muģ������ƽ������(�ɵ�),Pm����ʽ�������,muu����ʽ�������(�ɵ�),sigma��˹����ģ�͵ı�׼��(�ɵ�),
    %         selection_process�ɹ�ѡ��elitist��roulette wheel��Tournament,
    %         name����FCMʱ��ѡ���ת�������У�SFCM��WFCM��HTFCM��TLFCM
    %         options����matlab�����Ż�������ţ�ٷ���->����Ԥѧϰ�Ż���
    %         beta��data_error������probswap��λ��������
    % Output: data_MFEA������ʱ�䡢ÿ�����ֵ����ø�������������ʹ�õģ���Ϊÿ�����䣬�����������ûɶ���壩�������۴�����
    global Nk N gen Nc sdata g batch
    tic
    Ncc = Nc/batch;
    NN=0;
    %% 0.˥��ϵ���������Ÿ������
    alpha = zeros(Nk,Nc);
    Best = zeros(Nc,Nc,Nk);
    EvBestFitness = zeros(gen+1,batch,Nk);           %ÿ����õ���Ӧ��ֵ
    TotalEvaluations=zeros(gen+1,1,Nk);              %ÿ��ÿ���������۴���
    allEvBestFitness = zeros(gen+1,Nc,Nk);           %��������ÿ����õ���Ӧ��ֵ
    %% 1.������ѧϰNk��FCM
    for j = 1:Nk
        % 1.1 ��ʼ��FCM�����񲢼���alpha
        disp(['Nk = ',num2str(j)]);
        alpha = calalpha(alpha,sdata,Best(:,:,1:j-1),g);
        for bb = 1:Ncc       
            disp(['batch = ',num2str(bb)]);
            Task = TASK();
            Task = initTASK(Task,Nc,sdata,g,alpha(1:j,:),Best(:,:,1:j-1),batch,bb);

            % 1.2 ��ʼ����Ⱥ
            Population = INDIVIDUAL();                    %���ɳ�ʼ��Ⱥ
            Population = initPOP(Population,N,Task.D_multitask,Task.M);

            % 1.3 ���ݶ����񻷾��е�ÿ���Ż���������ÿ����������Ӵ���
            [Population,TotalEvaluations(1,1,j)] = evaluate(Population,Task);

            % 1.4 �����ʼ����Ⱥ�����صȼ��Լ���������
            [Population,EvBestFitness(1,:,j),bestind] = Calfactor(Population,beta);

            % 1.5 �Ż�����
            for i = 1:gen
                %4.1 ������콻��
                Offspring  = GA_MFEA(Population,rmp,Pc,Pm,Pd,mu,muu,sigma,probswap,i,gen);
                %4.2 �������Ӵ���
                [Offspring,TotalEvaluations(i+1,1,j)] = evaluate(Offspring,Task);
                TotalEvaluations(i+1,1,j) = TotalEvaluations(i+1,1,j) + TotalEvaluations(i+1,1,j);
                %4.3 ��Ⱥ�ϲ�
                intpopulation = combpop(Population,Offspring);
                %4.4 ���±�����Ӧ�ȣ��������أ����صȼ�
                [intpopulation,EvBestFitness(i+1,:,j),bestind] = Calfactor(intpopulation,beta);
                %4.5 ����ѡ��
                Population = EnvironmentalSelection(intpopulation,selection_process,N,Task.M);
                disp(['BMFEA Generation = ', num2str(i), ' EvBestFitness = ', num2str(EvBestFitness(i+1,:,j))]);%Ϊ�˼�¼��ʼ����ֵ���Դ���+1
            end

            % 1.6����ģ����֪ͼ�����۷�ʽ
            Best((1+batch*(bb-1)):(bb*batch),:,j) = bestind;
            allEvBestFitness(:,(1+batch*(bb-1)):(bb*batch),:) = EvBestFitness;
        end
        % 1.7���������˳�
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

    %% 2.��¼�㷨���
    data_BMFEA.Nk = NN;
    data_BMFEA.wall_clock_time=toc;
    data_BMFEA.alpha=alpha;
    data_BMFEA.EvBestFitness=EvBestFitness;
    data_BMFEA.allEvBestFitness=allEvBestFitness;
    data_BMFEA.bestInd_data=Best;
    data_BMFEA.TotalEvaluations=TotalEvaluations;
end