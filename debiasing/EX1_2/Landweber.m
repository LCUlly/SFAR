% solve equation by Landweber iteration method
% Ax=y;
% t^(1-theta) D^theta(t^(theta-1) x(t)) + A*A x(t) = A*y
function [x_land] = Landweber(y,N)
m = 100; tau = 1/m; s = (0:tau:1)'; t = 0:tau:1;
K = min(s,t).*(1-max(s,t)); 
x_real = (-6*t.^2.*(1-t).*(2-8*t+7*t.^2))';
A =tau*K*diag([1/2,ones(1,m-1),1/2]);

x = zeros(m+1,1);
g = A'*(y-A*x);
omg = 0.5/(norm(A,2))^2;

% N=10000;
Err_land=zeros(N+1,1); Err_land(1) = norm(x-x_real,2)/norm(x_real,2);% save error
Residual =zeros(N+1,1);
X = zeros(m+1,N+1);  % save x(t) in X
for n = 1:N
    x_old = x;
    x = x_old +omg*g; % update solution
    % x(1) = x(2);x(end) =x(end-1);
    X(:,n+1) = x;
    Err_land(n+1) = norm(x-x_real,2)/norm(x_real,2);
    g = A'*(y-A*x); % update search direction;

    Residual(n) = norm(A*x - y);

    x_land = x;


    if n>10 && Residual(n)>Residual(n-1)
        fprintf('dont decrease')
        % break;
    end
    % if Err_land(n+1)<min(Err_land(1:n))
    %     x_land = x;
    %     num_land = n;
    % end
    figure(1);
    subplot(2,1,1),plot(t',X(:,n+1),t,x_real,'r',LineWidth=3)
    % 创建 ylabel
    ylabel({'x'});
    % 创建 xlabel
    xlabel({'n'});
    % 创建 title
    title({'debias Landweber'});
    subplot(2,1,2),plot(max(1,n-100):n+1,Err_land(max(1,n-100):n+1))
    % 创建 title
    title({'L^2 error'});
end
