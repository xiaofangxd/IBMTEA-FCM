function obj = WFCM(var,a,b)
    %Waveletת������
    obj = (1-((var-a)/b).^2).*exp(-((var-a).^2)/(2*b^2));
end

