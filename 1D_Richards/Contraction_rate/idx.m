function index=idx(s,b_values)
%%% This function computes the array indices corresponding to the string s

    global du
    global u_start
    
    % s is assumed
    s=max(s,0);
    index = ceil((s-u_start)/du);
    index=max(index,1);
    index = min(index, length(b_values)-1);

end