function Fs_derivative = RdFds(b_iter,dbds_iter)

global lambda
[n, dim] = size(b_iter); 

  b_iter = max(b_iter,1e-70);

   T = (1-b_iter.^(1/lambda)).^lambda ;
   Y = (1-b_iter.^(1/lambda)).^(lambda-1);

Fs_derivative = zeros(size(b_iter));

   for i = 1 : n
   for j = 1: dim

    Fs_derivative(i,j) = 2.* (1-T(i,j)).*dbds_iter(i,j).*(b_iter(i,j))^((-lambda+1)/lambda)*sqrt(b_iter(i,j))*Y(i,j)+ ((dbds_iter(i,j)*(1-T(i,j))^2)/(2*sqrt(b_iter(i,j))));

    end
   end
end

