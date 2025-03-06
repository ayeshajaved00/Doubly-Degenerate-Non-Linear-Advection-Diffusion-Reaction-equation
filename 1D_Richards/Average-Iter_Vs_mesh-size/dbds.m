function bs_derivative = dbds(s,index_s,b_values)
    % Compute the derivative of the bs values with respect to s
    global du
    ispos = .5*(sign(s) + 1);    

    bs_derivative = ispos.*(b_values(index_s+1) - b_values(index_s)) / du;
end
