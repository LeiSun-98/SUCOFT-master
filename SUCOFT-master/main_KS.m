clc;
clear all;
close all;

addpath src;

load('data/simulated_correspondences_KS.mat');

outlier_ratios=[0.2 0.4 0.6 0.8 0.9 0.95 0.96 0.97 0.98 0.99];

for itr_outlier=2:size(store_R,2)

    outlier_ratio=outlier_ratios(itr_outlier-1);

    for nr_run=1 %1:size(store_R,1)

        % Getting correspondences
        R_gt=cell2mat(store_R(nr_run,itr_outlier));
        t_gt=cell2mat(store_t(nr_run,itr_outlier));
        pts_3d=cell2mat(store_pts_3d(nr_run,itr_outlier));
        pts_3d_=cell2mat(store_pts_3d_(nr_run,itr_outlier));

        % Set mininum inlier ratio
        min_inlier_ratio=0.01;  % at most 99% outliers

        % Warm-up
        [~, ~, ~] = SUCOFT_KS(pts_3d,pts_3d_,noise,min_inlier_ratio);

        % Running SUCOFT for known-scale PCR
        tic;
        [R_opt, t_opt, ~] = SUCOFT_KS(pts_3d,pts_3d_,noise,min_inlier_ratio);
        time(nr_run,itr_outlier)=toc();
        R_error(nr_run,itr_outlier)=AngularError(R_gt,R_opt)*180/pi;
        t_error(nr_run,itr_outlier)=norm(t_opt - t_gt');


        % Displaying the metrics
        disp(['Number of correspondences: ', num2str(n_ele)]);
        disp(['Noise: ', num2str(noise)]);
        disp(['Outlier ratio in percentage: ', num2str(outlier_ratio)]);
        disp(['Rotation error of SUCOFT in degrees: ', num2str(R_error(nr_run,itr_outlier))]);
        disp(['Translation error of SUCOFT in meters: ', num2str(t_error(nr_run,itr_outlier))]);
        disp(['Runtime of SUCOFT in seconds: ', num2str(time(nr_run,itr_outlier))]);
        disp('  ');


        % Plotting the correspondences
        figure;
        pc1=pointCloud(pts_3d(1:n_ele,:));
        pc2=pointCloud(pts_3d_(1:n_ele,:));
        pcshow([pc1.Location(:,1),pc1.Location(:,2),pc1.Location(:,3)],[0 0 1],'MarkerSize',70);
        hold on;
        pcshow([pc2.Location(:,1),pc2.Location(:,2),pc2.Location(:,3)],[0 1 1],'MarkerSize',70);
        for i=round(n_ele*outlier_ratio)+1:n_ele
            plot3([pts_3d(i,1),pts_3d_(i,1)],[pts_3d(i,2),pts_3d_(i,2)],[pts_3d(i,3),pts_3d_(i,3)],'g','LineWidth',3);
        end
        for i=1:round(n_ele*outlier_ratio)
            pp=plot3([pts_3d(i,1),pts_3d_(i,1)],[pts_3d(i,2),pts_3d_(i,2)],[pts_3d(i,3),pts_3d_(i,3)],'r','LineWidth',1);
            pp.Color(4) = 0.3;
        end
        grid off;
        axis off;

    end

end



