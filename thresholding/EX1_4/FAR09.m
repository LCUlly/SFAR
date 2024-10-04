% consider equation
% t^(1-theta) D^theta(t^(theta-1) x(t)) + A*A x(t) = A*y
function [x_threshold_09,x_09,num_09,Err_09,Err_threshold_09] = FAR09(noise_level)
thet = 0.9;  % thet in (0,1)
m = 100; tau = 1/m; s = (0:tau:1)'; t = 0:tau:1;
K = min(s,t).*(1-max(s,t));
t1 = (0:tau:1/5)'; t2 = (1/5+tau:tau:2/5)';
t3 = (2/5+tau:tau:3/5)';
t4 = (3/5+tau:tau:4/5)';t5 = (4/5+tau:tau:1)';
x_real = [zeros(size(t1));-ones(size(t2));zeros(size(t3));ones(size(t4));zeros(size(t5));];
A = tau*K*diag([1/2,ones(1,m-1),1/2]);% Rectangle integral formula
y_real = A*x_real;
% noise = (2*rand(size(y))-1);save('noise','noise');
load("noise.mat");%noise_level = 0.001;
y=y_real.*(1+noise*noise_level);% parameters in the equation


%  a concrete iterative regularization method.
h = 0.5/(norm(A,2))^2; % step size of artifitial time
N = 10000; % the maximum iterative number
% T = h*N; % maximum time to find
a = (1:N+1).^(1-thet)-(0:N).^(1-thet);

Residual =zeros(N+1,1);num=0;
X = zeros(m+1,N+1);  % save x(t) in X
Err_09=zeros(N+1,1); Err_09(1) = norm(X(:,1)-x_real,2)/norm(x_real,2);% save error
% Bx = f;
B = a(1)/gamma(2-thet)/h^thet*eye(m+1)...
    +1/2*(A'*A);
for n = 1:N
    f1 = 1/gamma(2-thet)/h^thet*X(:,2:n)*diag((n./(1:n-1)).^(1-thet))*(a(n-1:-1:1)-a(n:-1:2))'+...
        a(n)*X(:,1);
    f3 = -1/2*(A'*A)*X(:,n);
    f4 = A'*y;

    f = f1+f3+f4;
    x = B\f;
    X(:,n+1) = x;
    Err_09(n+1) = norm(x-x_real,2)/norm(x_real,2);

    Residual(n) = norm(A*x - y);
    % if n>=5 && Residual(n)/norm(y)<noise_level/2
    if n==N || Residual(n)/norm(y_real - y)<1.01
        num = num+1;
        if num ==1
            x_09 = x;
            num_09 = n;
            x_threshold_09 = x.*(1-(abs(x)<1e-1));
            Err_threshold_09 = norm(x_threshold_09-x_real,2)/norm(x_real,2);
            figure(3);
            %     subplot(2,1,1),plot(t',x,t,x_real,'r',LineWidth = 3)
            subplot(3,1,1),plot(t',x_threshold_09,t,x_real,'r')
            % 创建 ylabel
            ylabel({'x'});
            % 创建 xlabel
            xlabel({'n'});
            % 创建 title
            title({'\theta = 0.9'});
            subplot(3,1,2),plot(max(1,n-100):n+1,Err_09(max(1,n-100):n+1))
            % 创建 title
            title({'L^2 error'});
            subplot(3,1,3),plot(max(1,n-100):n,Residual(max(1,n-100):n)/norm(y_real - y))
            % 创建 title
            title({'\tau'});
            % pause
        end
    end
    if num>0 && n == num_09+50
        break;
    end

    if n>10 && Residual(n)>Residual(n-1)
        fprintf('dont decrease')
        % break;
    end

    figure(1);
    %     subplot(2,1,1),plot(t',X(:,n+1),t,x_real,'r',LineWidth=3)
    subplot(3,1,1),plot(t',x,t,x_real,'r')
    % 创建 ylabel
    ylabel({'x'});
    % 创建 xlabel
    xlabel({'n'});
    % 创建 title
    title({'\theta = 0.9'});
    subplot(3,1,2),plot(max(1,n-100):n+1,Err_09(max(1,n-100):n+1))
    % 创建 title
    title({'L^2 error'});
    subplot(3,1,3),plot(max(1,n-100):n,Residual(max(1,n-100):n)/norm(y_real - y))
    % 创建 title
    title({'\tau'});
end
