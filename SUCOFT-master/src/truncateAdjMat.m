function adj_in = truncateAdjMat(adj_in, thres_core_size)

%     adj_in = or((adj_in.*(adj_in*adj_in))>=thres_cliq_size, eye(size(adj_in,1)));

    adj_in = adj_in.*((adj_in*adj_in)>=(thres_core_size-1));
    
end
