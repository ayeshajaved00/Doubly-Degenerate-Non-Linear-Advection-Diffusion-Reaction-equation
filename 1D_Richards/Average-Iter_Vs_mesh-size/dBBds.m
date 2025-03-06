function Bs_derivative = dBBds(s,index_s,B_values)
    % Compute the derivative of the Bs values with respect to s
    global du
    ispos = .5*(sign(s) + 1);  
    Bs_derivative = 1- ispos + ispos.*(B_values(index_s+1) - B_values(index_s)) / du;
end
