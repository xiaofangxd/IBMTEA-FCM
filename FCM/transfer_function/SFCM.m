function obj = SFCM(var,g)
    %sigmoid×ª»»º¯Êý
    obj = 1./(1+exp(-g.*var));
end

