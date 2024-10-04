clc,clear
load("data_EX12.mat")

m = 100; tau = 1/m; s = (0:tau:1)'; t = 0:tau:1;
K = min(s,t).*(1-max(s,t));
A = tau*K*diag([1/2,ones(1,m-1),1/2]);
x_real = (-6*t.^2.*(1-t).*(2-8*t+7*t.^2))';
% 
%%  noise_level = 0.05;
[de_x_199_4_part] = FAR199(A*x_199_4,num_199_4);  
de_x_199_4 = 2*x_199_4 - de_x_199_4_part;%
de_err_199_4 = norm(de_x_199_4-x_real,2)/norm(x_real,2);%

[de_x_15_4_part] = FAR15(A*x_15_4,num_15_4);
de_x_15_4 = 2*x_15_4 - de_x_15_4_part;%
de_err_15_4 = norm(de_x_15_4-x_real,2)/norm(x_real,2);%

[de_x_land_4_part] = Landweber(A*x_land_4,num_land_4); 
de_x_land_4 = 2*x_land_4 - de_x_land_4_part;%
de_err_land_4 = norm(de_x_land_4-x_real,2)/norm(x_real,2);%

[de_x_09_4_part] = FAR09(A*x_09_4,num_09_4); 
de_x_09_4 = 2*x_09_4 - de_x_09_4_part;%
de_err_09_4 = norm(de_x_09_4-x_real,2)/norm(x_real,2);%

%% 
% noise_level = 0.01;  
[de_x_199_3_part] = FAR199(A*x_199_3,num_199_3);
de_x_199_3 = 2*x_199_3 - de_x_199_3_part;%
de_err_199_3 = norm(de_x_199_3-x_real,2)/norm(x_real,2);%

[de_x_15_3_part] = FAR15(A*x_15_3,num_15_3);
de_x_15_3 = 2*x_15_3 - de_x_15_3_part;%
de_err_15_3 = norm(de_x_15_3-x_real,2)/norm(x_real,2);%

[de_x_land_3_part] = Landweber(A*x_land_3,num_land_3); 
de_x_land_3 = 2*x_land_3 - de_x_land_3_part;%
de_err_land_3 = norm(de_x_land_3-x_real,2)/norm(x_real,2);%

[de_x_09_3_part] = FAR09(A*x_09_3,num_09_3); 
de_x_09_3 = 2*x_09_3 - de_x_09_3_part;%
de_err_09_3 = norm(de_x_09_3-x_real,2)/norm(x_real,2);%
%% 
% noise_level = 0.001;
[de_x_199_2_part] = FAR199(A*x_199_2,num_199_2);
de_x_199_2 = 2*x_199_2 - de_x_199_2_part;%
de_err_199_2 = norm(de_x_199_2-x_real,2)/norm(x_real,2);%

[de_x_15_2_part] = FAR15(A*x_15_2,num_15_2);
de_x_15_2 = 2*x_15_2 - de_x_15_2_part;%
de_err_15_2 = norm(de_x_15_2-x_real,2)/norm(x_real,2);%

[de_x_land_2_part] = Landweber(A*x_land_2,num_land_2); 
de_x_land_2 = 2*x_land_2 - de_x_land_2_part;%
de_err_land_2 = norm(de_x_land_2-x_real,2)/norm(x_real,2);%

[de_x_09_2_part] = FAR09(A*x_09_2,num_09_2); 
de_x_09_2 = 2*x_09_2 - de_x_09_2_part;%
de_err_09_2 = norm(de_x_09_2-x_real,2)/norm(x_real,2);%
%%
% noise_level = 0.0001;
[de_x_199_1_part] = FAR199(A*x_199_1,num_199_1);
de_x_199_1 = 2*x_199_1 - de_x_199_1_part;%
de_err_199_1 = norm(de_x_199_1-x_real,2)/norm(x_real,2);%

[de_x_15_1_part] = FAR15(A*x_15_1,num_15_1);
de_x_15_1 = 2*x_15_1 - de_x_15_1_part;%
de_err_15_1 = norm(de_x_15_1-x_real,2)/norm(x_real,2);%

[de_x_land_1_part] = Landweber(A*x_land_1,num_land_1); 
de_x_land_1 = 2*x_land_1 - de_x_land_1_part;%
de_err_land_1 = norm(de_x_land_1-x_real,2)/norm(x_real,2);%

[de_x_09_1_part] = FAR09(A*x_09_1,num_09_1); 
de_x_09_1 = 2*x_09_1 - de_x_09_1_part;%
de_err_09_1 = norm(de_x_09_1-x_real,2)/norm(x_real,2);%

% 
save('debiased_data_EX12')
% 
% 
% 
% 
