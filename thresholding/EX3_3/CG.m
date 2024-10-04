% solve equation by CG method
% refer to book <computational methods for inverse problems> page 32
function [x_threshold_CG,x_CG,num_CG,Err_CG,Err_threshold_CG] = CG(noise_level)
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
A = tau1*tau2*K;% Rectangle integral formula
xx_real = ((t1'-1/2).^2+(t2-1/2).^2<=9/64).*(4*sin(pi*t1').*sin(pi*t2)); 
[T1,T2] = meshgrid(t1',t2);
figure(1); mesh(T2,T1,xx_real);
x_real =  reshape(xx_real,(m1+1)*(m2+1),1);
y_real = A*x_real; y_real = reshape(y_real,(n1+1)*(n2+1),1);

% noise = (2*rand(size(y_real))-1);save('noise','noise');
load("noise.mat");
% noise_level = 0.001;
y=y_real.*(1+noise*noise_level);% parameters in the equation

N=1000;
X = zeros((m1+1)*(m2+1),N+1);   % save x(t) in X
x = X(:,1);
g = A'*(A*x-y);
p = -g;
dlt = norm(g,2)^2;

Err_CG=zeros(N+1,1); Err_CG(1) = norm(x-x_real,2)/norm(x_real,2);% save error
Residual =zeros(N+1,1);num=0;
for n = 1:N
    h = A'*A*p; tau = dlt/(p'*h);
    x_old = x;
    x = x_old +tau*p; % update solution
    X(:,n+1) = x;
    Err_CG(n+1) = norm(x-x_real,2)/norm(x_real,2);
    g = g+tau*h; % update gradient;
    dlt_old = dlt;
    dlt = norm(g,2)^2; beta = dlt/dlt_old;
    p = -g+beta*p; % update search direction;

    Residual(n) = norm(A*x - y);
    if  Residual(n)/norm(y_real - y)<1.01
        num = num+1;
        if num ==1
            x_CG = x;
            num_CG = n;
            x_threshold_CG = x.*(1-(abs(x)<1));
            Err_threshold_CG = norm(x_threshold_CG-x_real,2)/norm(x_real,2);
            xx = reshape(x_threshold_CG,m1+1,m2+1);
            [T1,T2] = meshgrid(t1',t2);
            figure(1); mesh(T2,T1,xx);
            % 创建 ylabel
            ylabel({'y'});
            % 创建 xlabel
            xlabel({'x'});
            % 创建 zlabel
            zlabel({'f(x,y)'});
            % 创建 title
            title({'CG'});
            figure(3); subplot(2,1,1),plot(max(1,n-100):n+1,Err_CG(max(1,n-100):n+1))
            subplot(2,1,2),plot(max(1,n-100):n,Residual(max(1,n-100):n)/norm(y_real - y))
            title({'L^2 error'});
            % pause
        end
    end
    if num>0 && n == num_CG+50
        break;
    end
    xx = reshape(x,m1+1,m2+1);
    [T1,T2] = meshgrid(t1',t2);
    figure(2); mesh(T2,T1,xx);
    % 创建 ylabel
    ylabel({'y'});
    % 创建 xlabel
    xlabel({'x'});
    % 创建 zlabel
    zlabel({'f(x,y)'});
    % 创建 title
    title({'CG'});
    figure(3); subplot(2,1,1),plot(max(1,n-100):n+1,Err_CG(max(1,n-100):n+1))
    subplot(2,1,2),plot(max(1,n-100):n,Residual(max(1,n-100):n)/norm(y_real - y))
    title({'L^2 error'});
end
