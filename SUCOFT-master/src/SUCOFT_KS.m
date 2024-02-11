function [R_opt, t_opt, best_set] = algorithm_KS(pts_3d,pts_3d_,noise,min_inlier_ratio)

bound=2*3.5*noise;

n_ele=size(pts_3d,1);

% adj_mat=diag(ones(1,n_ele));
adj_mat=zeros(n_ele,n_ele);
for i=1:n_ele-1
    v1=pts_3d(i,:)-pts_3d(i+1:end,:);
    v2=pts_3d_(i,:)-pts_3d_(i+1:end,:);
    adj_mat(i,i+1:n_ele)=(abs(sqrt(sum(v1.^2,2))-sqrt(sum(v2.^2,2)))'<=bound);
end

adj_mat=adj_mat+adj_mat';  %-diag(ones(1,n_ele));


tt=tic;

max_trunc_itr=100;

thres_core_size=max([5, round(min_inlier_ratio*n_ele)])-1;

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

[R_opt, t_opt, best_set] = robustSolver_KS(adj_, pts_3d, pts_3d_, noise);

end
