function obj = HTFCM(var,lambda)
    %Hyperbolic tangentת������
    obj = (exp(lambda.*var)-exp(-lambda.*var))./(exp(lambda.*var)+exp(-lambda.*var));
end

