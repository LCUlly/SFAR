close all
load("data_EX33_compare_CG.mat")
figure(3);
semilogx(...
    1:min(10001,num_land_3+51),Err_land_3(1:min(10001,num_land_3+51)),'-',...
    1:min(10001,num_09_3+51),Err_09_3(1:min(10001,num_09_3+51)),'-.',...
    1:num_15_3+51,Err_15_3(1:num_15_3+51),'--', ...
    1:num_199_3+51,Err_199_3(1:num_199_3+51),'o-', ...
    1:num_CG_3+51,Err_CG_3(1:num_CG_3+51),'-', ...
    LineWidth=2)
% 创建 title
title({'L^2 errors for Example 3(3)'},fontsize=20);
legend(...%'CG', ...
    'landweber','FAR, \theta = 0.9',...
    'FAR, \theta = 1.5', ... %'FAR, \theta = 1.9', ...
    'FAR, \theta = 1.99', ...
    'CG',...
    fontsize=20)

m1 = 50; tau1 = 1/m1; % m1=n1
t1 = 0:tau1:1;
m2 = 50; tau2 = 1/m2; % n2 = m2;
t2 = 0:tau2:1;
xx_real = ((t1'-1/2).^2+(t2-1/2).^2<=9/64).*(4*sin(pi*t1').*sin(pi*t2)); %good
[T1,T2] = meshgrid(t1',t2);
figure(1);
subplot(2,3,1),mesh(T2,T1,xx_real);
% 创建 ylabel
ylabel({'y'});
% 创建 xlabel
xlabel({'x'});
% 创建 zlabel
zlabel({'f(x,y)'});
% 创建 title
title({'exact f(x,y)'});

subplot(2,3,2),
xx = reshape(x_land_3,m1+1,m2+1);
% contour(T2,T1,xx-xx_real,'ShowText','on');
mesh(T2,T1,xx);
% 创建 ylabel
ylabel({'y'});
% 创建 xlabel
xlabel({'x'});
% 创建 zlabel
zlabel({'f(x,y)'});
% 创建 title
title({'Landweber'});

subplot(2,3,3),
xx = reshape(x_CG_3,m1+1,m2+1);
% contour(T2,T1,xx-xx_real,'ShowText','on');
mesh(T2,T1,xx);
% 创建 ylabel
ylabel({'y'});
% 创建 xlabel
xlabel({'x'});
% 创建 zlabel
zlabel({'f(x,y)'});
% 创建 title
title({'CG'});

subplot(2,3,4),
xx = reshape(x_threshold_land_3,m1+1,m2+1);
% contour(T2,T1,xx-xx_real,'ShowText','on');
mesh(T2,T1,xx);
% 创建 ylabel
ylabel({'y'});
% 创建 xlabel
xlabel({'x'});
% 创建 zlabel
zlabel({'f(x,y)'});
% 创建 title
title({'thresholded Landweber (\theta = 1)'});

subplot(2,3,5),
xx = reshape(x_threshold_15_3,m1+1,m2+1);
% contour(T2,T1,xx-xx_real,'ShowText','on');
mesh(T2,T1,xx);
% 创建 ylabel
ylabel({'y'});
% 创建 xlabel
xlabel({'x'});
% 创建 zlabel
zlabel({'f(x,y)'});
% 创建 title
title({'\theta = 1.5'});

subplot(2,3,6),
xx = reshape(x_threshold_199_3,m1+1,m2+1);
% contour(T2,T1,xx-xx_real,'ShowText','on');
mesh(T2,T1,xx);
% mesh(T2,T1,xx);
% 创建 ylabel
ylabel({'y'});
% 创建 xlabel
xlabel({'x'});
% 创建 zlabel
zlabel({'f(x,y)'});
% 创建 title
title({'\theta = 1.99'});

