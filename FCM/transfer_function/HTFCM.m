function obj = HTFCM(var,lambda)
    %Hyperbolic tangent×ª»»º¯Êý
    obj = (exp(lambda.*var)-exp(-lambda.*var))./(exp(lambda.*var)+exp(-lambda.*var));
end

