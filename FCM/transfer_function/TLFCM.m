function obj = TLFCM(var,t)
    %Threshold linearת������
    obj = var;
    obj(obj<t | obj==t) = 0;
    obj(obj>t) = obj(obj>t) - t;
end

