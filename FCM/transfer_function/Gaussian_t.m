function obj = Gaussian_t(var,TT,DD)
    %Gaussianת������
    obj = exp(-(var-TT).^2./(2*DD.^2));
end

