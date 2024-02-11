function adj_ = iterativeAdjMatRefining(adj_,thres_core_size,max_trunc_itr)

for rep=1:max_trunc_itr
    
    adj_last=adj_;
    
    adj_=truncateAdjMat(adj_, thres_core_size);
    
    if max(max(abs(adj_-adj_last)))==0
        return
    end
    
end

end


