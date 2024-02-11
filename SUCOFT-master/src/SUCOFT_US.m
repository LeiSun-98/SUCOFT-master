function [scale_opt, R_opt, t_opt, best_set] = algorithm_US(pts_3d,pts_3d_,noise,min_inlier_ratio)

bound=4*noise;

max_trunc_itr=100;

n_ele=size(pts_3d,1);

scale_mat=zeros(n_ele,n_ele);
length_mat=zeros(n_ele,n_ele);
for i=1:n_ele
    v1=pts_3d(i,:)-pts_3d(i+1:n_ele,:);
    v2=pts_3d_(i,:)-pts_3d_(i+1:n_ele,:);
    length_mat(i,i+1:n_ele)=1./sqrt(sum(v1.^2,2))';
    scale_mat(i,i+1:n_ele)=sqrt(sum(v2.^2,2))'.*length_mat(i,i+1:n_ele);
end

scale_mat=scale_mat+scale_mat';
length_mat=length_mat+length_mat';
length_mat=bound*length_mat;

tt=tic;

thres_core_size=max([5, round(min_inlier_ratio*n_ele)])-1;

adj_mat=buildAdjMat_US(n_ele, thres_core_size, scale_mat, length_mat);

adj_ = iterativeAdjMatRefining(adj_mat,thres_core_size,max_trunc_itr);

adj_mat_refined=adj_;

if sum(sum(adj_mat_refined))~=thres_core_size*(thres_core_size+1)

    check_vec=sort(sum(adj_mat_refined),'descend');
    for i=1:n_ele
        if check_vec(i)+1<i
            cliq_size_candidates=i-2:-1:thres_core_size;
            break
        end
    end

    max(cliq_size_candidates)
    
    for poss_cliq_size=cliq_size_candidates

        thres_core_size=poss_cliq_size;
        
        adj_ = iterativeAdjMatRefining(adj_mat_refined,thres_core_size,max_trunc_itr);

        if sum(sum(adj_))>0
            break
        end
        
    end
    
end

disp(['K-super-core size is: ',num2str(thres_core_size)]);

toc(tt);

[scale_opt, R_opt, t_opt, best_set] = robustSolver_US(adj_, pts_3d, pts_3d_, noise);

end


