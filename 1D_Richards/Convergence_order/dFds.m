function Fs_derivative = dFds(b_iter,dbds_iter)

    global lambda

  b_iter = max(b_iter,1e-9);

   T = (1-b_iter.^(1/lambda)).^lambda ;
   Y = (1-b_iter.^(1/lambda)).^(lambda-1);

Fs_derivative = zeros(1,length(b_iter)); 

   for i = 1 : length(b_iter)

    Fs_derivative(i) = 2.* (1-T(i)).*dbds_iter(i).*(b_iter(i))^((-lambda+1)/lambda)*sqrt(b_iter(i))*Y(i)+ ((dbds_iter(i)*(1-T(i))^2)/(2*sqrt(b_iter(i))));

    end
end

