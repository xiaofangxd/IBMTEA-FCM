function obj = Gaussian_t(var,TT,DD)
    %Gaussian×ª»»º¯Êý
    obj = exp(-(var-TT).^2./(2*DD.^2));
end

