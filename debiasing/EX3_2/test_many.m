clc,clear
load("data_EX32.mat")

n1 = 50; sau1 = 1/n1; m1 = 50; tau1 = 1/m1; % m1=n1
s1 = (0:sau1:1)'; t1 = 0:tau1:1;
n2 = 50; sau2 = 1/n2;m2 = 50; tau2 = 1/m2; % n2 = m2;
s2 = (0:sau2:1)'; t2 = 0:tau2:1;
K = zeros((n1+1)*(n2+1),(m1+1)*(m2+1));

for i = 1:n1+1
    for j = 1:n2+1
        KK =  sin(pi*(s1(i)-t1').^2).*sin(pi*(s2(j)-t2).^2);
        KK = diag([1/2,ones(1,m1-1),1/2])*KK*diag([1/2,ones(1,m2-1),1/2]);
       K(i+(j-1)*(n2+1),:) = reshape(KK,1,(m1+1)*(m2+1)); 
    end
end
A = tau1*tau2*K;%  integral formula% 
xx_real =  1+ (t1').^2.*t2.^2;
x_real =  reshape(xx_real,(m1+1)*(m2+1),1);
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
% % noise_level = 0.001;
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
% % noise_level = 0.0001;
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
save('debiased_data_EX32')
% 
% 
% 
% 
