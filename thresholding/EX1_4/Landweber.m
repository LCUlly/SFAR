% solve equation by Landweber iteration method
% Ax=y;
% t^(1-theta) D^theta(t^(theta-1) x(t)) + A*A x(t) = A*y
function [x_threshold_land,x_land,num_land,Err_land,Err_threshold_land] = Landweber(noise_level)
m = 100; tau = 1/m; s = (0:tau:1)'; t = 0:tau:1;
K = min(s,t).*(1-max(s,t));
t1 = (0:tau:1/5)'; t2 = (1/5+tau:tau:2/5)';
t3 = (2/5+tau:tau:3/5)';
t4 = (3/5+tau:tau:4/5)';t5 = (4/5+tau:tau:1)';
x_real = [zeros(size(t1));-ones(size(t2));zeros(size(t3));ones(size(t4));zeros(size(t5));];
A = tau*K*diag([1/2,ones(1,m-1),1/2]);% Rectangle integral formula
y_real = A*x_real;
% noise = (2*rand(size(y_real))-1);save('noise','noise');
load("noise.mat");%noise_level = 0.001;
y=y_real.*(1+noise*noise_level);% parameters in the equation


x = zeros(m+1,1);
g = A'*(y-A*x);
omg = 0.5/(norm(A,2))^2;

N=10000;Err_land=zeros(N+1,1); Err_land(1) = norm(x-x_real,2)/norm(x_real,2);% save error
Residual =zeros(N+1,1);num=0;
X = zeros(m+1,N+1);  % save x(t) in X
for n = 1:N
    x_old = x;
    x = x_old +omg*g; % update solution
    % x(1) = x(2);x(end) =x(end-1);
    X(:,n+1) = x;
    Err_land(n+1) = norm(x-x_real,2)/norm(x_real,2);
    g = A'*(y-A*x); % update search direction;

    Residual(n) = norm(A*x - y);
    % if n>=5 && Residual(n)/norm(y)<noise_level/1
    if n==N || Residual(n)/norm(y_real - y)<1.01
        num = num+1;
        if num ==1
            x_land = x;
            num_land = n;
            x_threshold_land = x.*(1-(abs(x)<1e-1));
            Err_threshold_land = norm(x_threshold_land-x_real,2)/norm(x_real,2);
            figure(3);
            %     subplot(2,1,1),plot(t',x,t,x_real,'r',LineWidth = 3)
            subplot(3,1,1),plot(t',x_threshold_land,t,x_real,'r')
            % 创建 ylabel
            ylabel({'x'});
            % 创建 xlabel
            xlabel({'n'});
            % 创建 title
            title({'\theta = 0.9'});
            subplot(3,1,2),plot(max(1,n-100):n+1,Err_land(max(1,n-100):n+1))
            % 创建 title
            title({'L^2 error'});
            subplot(3,1,3),plot(max(1,n-100):n,Residual(max(1,n-100):n)/norm(y_real - y))
            % 创建 title
            title({'\tau'});
            % pause
        end
    end
    if num>0 && n == num_land+50
        break;
    end



    figure(1);
    %     subplot(2,1,1),plot(t',X(:,n+1),t,x_real,'r',LineWidth=3)
    subplot(3,1,1),plot(t',x,t,x_real,'r')
    % 创建 ylabel
    ylabel({'x'});
    % 创建 xlabel
    xlabel({'n'});
    % 创建 title
    title({'Landweber'});
    subplot(3,1,2),plot(max(1,n-100):n+1,Err_land(max(1,n-100):n+1))
    % 创建 title
    title({'L^2 error'});
    subplot(3,1,3),plot(max(1,n-100):n,Residual(max(1,n-100):n)/norm(y_real - y))
    % 创建 title
    title({'\tau'});
end
