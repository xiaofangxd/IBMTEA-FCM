function obj = SFCM(var,g)
    %sigmoidת������
    obj = 1./(1+exp(-g.*var));
end

