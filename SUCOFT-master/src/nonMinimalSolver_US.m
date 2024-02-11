function [s_best,R_opt,t_opt] = nonMinimalSolver_US(pts_3d,pts_3d_,best_set)

q_=zeros(3,1);
p_=zeros(3,1);

best_size=length(best_set);

if best_size<=1
    R_opt=eye(3);
    t_opt=zeros(3,1);
    s_best=1;
    return
end

for i=1:best_size

    q_=q_+pts_3d_(best_set(i),:)';

    p_=p_+pts_3d(best_set(i),:)';

end

p_=p_/best_size;
q_=q_/best_size;

s_=zeros(1,best_size);
s_scale=zeros(1,best_size);

for i=1:best_size
    s_(i)=norm(pts_3d_(best_set(i),:)'-q_)/norm(pts_3d(best_set(i),:)'-p_);
    s_scale(i)=norm(pts_3d(best_set(i),:)'-p_);
end

s_best=0;s_sum=0;

for i=1:best_size
s_best=s_best+s_(i)/(0.01/(s_scale(i))^2);
s_sum=s_sum+1/(0.01/s_scale(i)^2);
end

s_best=s_best/s_sum;

H=zeros(3,3);
for i=1:best_size
    H=H+(pts_3d(best_set(i),:)'-p_)*(pts_3d_(best_set(i),:)'-q_)';
end

[U,~,V]=svd(H);

R_opt=V*U';

t_opt=q_-s_best*R_opt*p_;


end

