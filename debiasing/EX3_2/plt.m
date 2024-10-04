
close all
load('debiased_data_EX32')

n1 = 50; sau1 = 1/n1; m1 = 50; tau1 = 1/m1; % m1=n1
s1 = (0:sau1:1)'; t1 = 0:tau1:1;
n2 = 50; sau2 = 1/n2;m2 = 50; tau2 = 1/m2; % n2 = m2;
s2 = (0:sau2:1)'; t2 = 0:tau2:1;
xx_real =  1+ (t1').^2.*t2.^2;
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
xx = reshape(de_x_land_3,m1+1,m2+1);
mesh(T2,T1,xx);
% 创建 ylabel
ylabel({'y'});
% 创建 xlabel
xlabel({'x'});
% 创建 zlabel
zlabel({'f(x,y)'});
% 创建 title
title({'debias Landweber'});

subplot(2,3,4),
xx = reshape(de_x_09_3,m1+1,m2+1);
mesh(T2,T1,xx);
% 创建 ylabel
ylabel({'y'});
% 创建 xlabel
xlabel({'x'});
% 创建 zlabel
zlabel({'f(x,y)'});
% 创建 title
title({'debias \theta = 0.9'});

subplot(2,3,5),
xx = reshape(de_x_15_3,m1+1,m2+1);
mesh(T2,T1,xx);
% 创建 ylabel
ylabel({'y'});
% 创建 xlabel
xlabel({'x'});
% 创建 zlabel
zlabel({'f(x,y)'});
% 创建 title
title({'debias \theta = 1.5'});

subplot(2,3,6),
xx = reshape(de_x_199_3,m1+1,m2+1);
mesh(T2,T1,xx);
% 创建 ylabel
ylabel({'y'});
% 创建 xlabel
xlabel({'x'});
% 创建 zlabel
zlabel({'f(x,y)'});
% 创建 title
title({'debias \theta = 1.99'});
