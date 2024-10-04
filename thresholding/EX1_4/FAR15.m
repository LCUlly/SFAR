% consider equation
% t^(1-theta) D^theta(t^(theta-1) x(t)) + A*A x(t) = A*y
function [x_threshold_15,x_15,num_15,Err_15,Err_threshold_15] = FAR15(noise_level)
thet = 1.5;  % thet in (0,1)
m = 100; tau = 1/m; s = (0:tau:1)'; t = 0:tau:1;
K = min(s,t).*(1-max(s,t));
t1 = (0:tau:1/5)'; t2 = (1/5+tau:tau:2/5)';
t3 = (2/5+tau:tau:3/5)';
t4 = (3/5+tau:tau:4/5)';t5 = (4/5+tau:tau:1)';
x_real = [zeros(size(t1));-ones(size(t2));zeros(size(t3));ones(size(t4));zeros(size(t5));];
A = tau*K*diag([1/2,ones(1,m-1),1/2]);% Rectangle integral formula
y_real = A*x_real;
% noise = (2*rand(size(y))-1); save('noise','noise');
load("noise.mat");%noise_level = 0.001;
y=y_real.*(1+noise*noise_level);% parameters in the equation


%  a concrete iterative regularization method.
h = 0.5/(norm(A,2))^2; % step size of artifitial time
N = 10000; % the maximum iterative number
% T = h*N; % maximum time to find
b = (1:N+1).^(2-thet)-(0:N).^(2-thet);

Residual =zeros(N+1,1);num=0;
X = zeros(m+1,N+1);  % save x(t) in X
Err_15=zeros(N+1,1); Err_15(1) = norm(X(:,1)-x_real,2)/norm(x_real,2);% save error
% Bx = f;
B = b(1)/gamma(3-thet)/h^thet*eye(m+1)...
    +1/2*(A'*A);
for n = 1:N
    f1 = b(1)/gamma(3-thet)/h^thet*((n-1)/n)^(thet-1)*X(:,n)...
        +1/gamma(3-thet)/h^thet*(X(:,2:n)*diag(((1:n-1)./n).^(thet-1))-X(:,1:n-1)*diag(((0:n-2)./n).^(thet-1)))*(b(n-1:-1:1)-b(n:-1:2))';

    f3 = -1/2*(A'*A)*X(:,n);
    f4 = A'*y;

    f = f1+f3+f4;
    x = B\f;
    X(:,n+1) = x;

    Err_15(n+1) = norm(x-x_real,2)/norm(x_real,2);
    Residual(n) = norm(A*x - y);
    % if n>=5 && Residual(n)/norm(y)<0.9*noise_level
    if n==N || Residual(n)/norm(y_real - y)<1.01
        num = num+1;
        if num ==1
            x_15 = x;
            num_15 = n;
            x_threshold_15 = x.*(1-(abs(x)<1e-1));
            Err_threshold_15 = norm(x_threshold_15-x_real,2)/norm(x_real,2);
            figure(3);
            %     subplot(2,1,1),plot(t',x,t,x_real,'r',LineWidth = 3)
            subplot(3,1,1),plot(t',x_threshold_15,t,x_real,'r')
            % 创建 ylabel
            ylabel({'x'});
            % 创建 xlabel
            xlabel({'n'});
            % 创建 title
            title({'\theta = 1.5'});
            subplot(3,1,2),plot(max(1,n-100):n+1,Err_15(max(1,n-100):n+1))
            % 创建 title
            title({'L^2 error'});
            subplot(3,1,3),plot(max(1,n-100):n,Residual(max(1,n-100):n)/norm(y_real - y))
            % 创建 title
            title({'\tau'});
            % pause
        end
    end
    if num>0 && n == num_15+50
        break;
    end

    if n>10 && Residual(n)>Residual(n-1)
        fprintf('dont decrease')
        % break;
    end

    figure(3);
    %     subplot(2,1,1),plot(t',x,t,x_real,'r',LineWidth=3)
    subplot(3,1,1),plot(t',x,t,x_real,'r')
    % 创建 ylabel
    ylabel({'x'});
    % 创建 xlabel
    xlabel({'n'});
    % 创建 title
    title({'\theta = 1.5'});
    subplot(3,1,2),plot(max(1,n-100):n+1,Err_15(max(1,n-100):n+1))
    % 创建 title
    title({'L^2 error'});
    subplot(3,1,3),plot(max(1,n-100):n,Residual(max(1,n-100):n)/norm(y_real - y))
    % 创建 title
    title({'\tau'});
end
