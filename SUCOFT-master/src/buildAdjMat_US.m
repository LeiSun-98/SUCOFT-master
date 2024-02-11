function adj_mat=buildAdjMat_US(n_ele, thres_core_size, scale_mat, length_mat)

adj_mat=zeros(n_ele,n_ele);
for i=1:n_ele-1
    for j=i+1:n_ele
        scale_this=scale_mat(i,j);
        len_this=length_mat(i,j);
        vec1=scale_mat(i,:);
        vec2=scale_mat(j,:);
        vec1(i)=scale_this;
        vec2(j)=scale_this;
        eligibility=0;
        for k=1:n_ele
            if abs(vec1(k)-scale_this)<=length_mat(i,k)+len_this && k~=i && k~=j
                if abs(vec2(k)-scale_this)<=length_mat(j,k)+len_this
                    if abs(vec1(k)-vec2(k))<=length_mat(i,k)+length_mat(j,k)
                        eligibility=eligibility+1;
                        if eligibility>=thres_core_size-1
                            break
                        end
                    end
                end
            end
        end
        if eligibility>=thres_core_size-1
            adj_mat(i,j)=1;
            adj_mat(j,i)=1;
        end
    end
end

end