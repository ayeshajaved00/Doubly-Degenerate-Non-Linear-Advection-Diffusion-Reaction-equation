function Bs = BB(s,index_s,B_values)
    % Perform linear interpolation to estimate the value of B for a given s

    global du
    global u_start

    sneg=min(s, 0);
    spos=max(s,0);
    Bs= sneg + (B_values(index_s).*(u_start+du*index_s-spos)+ B_values(index_s+1).*(spos-(u_start+du*(index_s-1))))/du;

end