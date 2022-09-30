function [objective,rnvec,calls] = CalObj(Task,rnvec,i,skill_factor)
% ������Ⱥ���и����ڵ�i�������Ŀ�꺯��ֵ��factorial_costs��
% Input:
% Task������Ϣ��ά�ȣ����½磬��������rnvec�����۵�Ⱦɫ�塢i���۵�i������skill_factor������Ⱦɫ��ļ�������
% Output:
% objective������Ⱦɫ���Ŀ�꺯��ֵ��rnvec����Ԥѧϰ���磬�޸ĺ��Ⱦɫ�塢calls���۴���
    global N
    calind = find(skill_factor == 0 | skill_factor == i);%�ҵ���Ҫ���۵���Ⱥ����
    objective = inf*ones(N,1);%Ŀ�꺯��ֵ��ʼ��
    d = Task.Tdims(i);
    nvars = rnvec(calind,1:d);
    x=nvars;
    objective1=Task.fun(i).fnc(x);
    calls = length(calind);%��Ӧ�����۴���
    objective(calind)=objective1;
end