function [R_opt,t_opt] = nonMinimalSolver_KS(pts_3d,pts_3d_,best_set)

q_=zeros(3,1);
p_=zeros(3,1);

best_size=length(best_set);

if best_size<=1
    R_opt=eye(3);
    t_opt=zeros(3,1);
    return
end

for i=1:best_size

    q_=q_+pts_3d_(best_set(i),:)';

    p_=p_+pts_3d(best_set(i),:)';

end

p_=p_/best_size;
q_=q_/best_size;

H=zeros(3,3);
for i=1:best_size
    H=H+(pts_3d(best_set(i),:)'-p_)*(pts_3d_(best_set(i),:)'-q_)';
end

[U,~,V]=svd(H);

R_opt=V*U';

t_opt=q_-R_opt*p_;


end

