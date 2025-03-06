function bs = bs(s,index_s,b_values)

    % Perform linear interpolation to estimate the value of b for a given s
    global du
    global u_start

    spos=max(s,0);
    bs = (b_values(index_s).*(u_start+du*index_s-spos)+ b_values(index_s+1).*(spos-(u_start+du*(index_s-1))))/du;
    bs=min(bs,1);
end
   