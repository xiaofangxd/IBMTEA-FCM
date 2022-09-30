function [objective,rnvec,calls] = CalObj(Task,rnvec,i,skill_factor)
% 计算种群所有个体在第i个任务的目标函数值（factorial_costs）
% Input:
% Task任务信息（维度，上下界，函数）、rnvec被评价的染色体、i评价第i个任务、skill_factor被评价染色体的技能因子
% Output:
% objective被评价染色体的目标函数值、rnvec由于预学习出界，修改后的染色体、calls评价次数
    global N
    calind = find(skill_factor == 0 | skill_factor == i);%找到需要评价的种群个体
    objective = inf*ones(N,1);%目标函数值初始化
    d = Task.Tdims(i);
    nvars = rnvec(calind,1:d);
    x=nvars;
    objective1=Task.fun(i).fnc(x);
    calls = length(calind);%适应度评价次数
    objective(calind)=objective1;
end