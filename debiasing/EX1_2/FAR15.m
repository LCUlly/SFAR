% consider equation
% t^(1-theta) D^theta(t^(theta-1) x(t)) + A*A x(t) = A*y
function [x_15] = FAR15(y,N)
thet = 1.5;  % thet in (0,1)
m = 100; tau = 1/m; s = (0:tau:1)'; t = 0:tau:1;
K = min(s,t).*(1-max(s,t));
x_real = (-6*t.^2.*(1-t).*(2-8*t+7*t.^2))';
A = tau*K*diag([1/2,ones(1,m-1),1/2]);

%  a concrete iterative regularization method.
h = 0.5/(norm(A,2))^2; % step size of artifitial time
% N = 10000; % the maximum iterative number
% T = h*N; % maximum time to find
b = (1:N+1).^(2-thet)-(0:N).^(2-thet);

Residual =zeros(N+1,1);
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


    x_15 = x;

    if n>10 && Residual(n)>Residual(n-1)
        fprintf('dont decrease')
        % break;
    end

    figure(1);
    subplot(2,1,1),plot(t',x,t,x_real,'r',LineWidth=3)
    % 创建 ylabel
    ylabel({'x'});
    % 创建 xlabel
    xlabel({'n'});
    % 创建 title
    title({'debias \theta = 1.5'});
    subplot(2,1,2),plot(max(1,n-100):n+1,Err_15(max(1,n-100):n+1))
    % 创建 title
    title({'L^2 error'});
end
