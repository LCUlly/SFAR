% consider equation
% t^(1-theta) D^theta(t^(theta-1) x(t)) + A*A x(t) = A*y
function [x_199] = FAR199(y,N)
thet = 1.99;  % thet in (0,1)
n1 = 50; sau1 = 1/n1; m1 = 50; tau1 = 1/m1; % m1=n1
s1 = (0:sau1:1)'; t1 = 0:tau1:1;
n2 = 50; sau2 = 1/n2;m2 = 50; tau2 = 1/m2; % n2 = m2;
s2 = (0:sau2:1)'; t2 = 0:tau2:1;
K = zeros((n1+1)*(n2+1),(m1+1)*(m2+1));
for i = 1:n1+1
    for j = 1:n2+1
        KK = sin(pi*(s1(i)-t1').^2).*sin(pi*(s2(j)-t2).^2);
        KK = diag([1/2,ones(1,m1-1),1/2])*KK*diag([1/2,ones(1,m2-1),1/2]);
       K(i+(j-1)*(n2+1),:) = reshape(KK,1,(m1+1)*(m2+1)); 
    end
end
A = tau1*tau2*K;% Rectangle integral formula
xx_real =  1+ (t1').^2.*t2.^2;
[T1,T2] = meshgrid(t1',t2);
figure(1); mesh(T2,T1,xx_real);
x_real =  reshape(xx_real,(m1+1)*(m2+1),1);

%  a concrete iterative regularization method.
h = 1/(norm(A,2))^2; % step size of artifitial time
% T = h*N; % maximum time to find
b = (1:N+1).^(2-thet)-(0:N).^(2-thet);

Residual =zeros(N+1,1);
X = zeros((m1+1)*(m2+1),N+1);  % save x(t) in X
Err_199=zeros(N+1,1); Err_199(1) = norm(X(:,1)-x_real,2)/norm(x_real,2);% save error
% Bx = f;
B = b(1)/gamma(3-thet)/h^thet*eye((m1+1)*(m2+1))...
    +1/2*(A'*A);
for n = 1:N
    f1 = b(1)/gamma(3-thet)/h^thet*((n-1)/n)^(thet-1)*X(:,n)...
        +1/gamma(3-thet)/h^thet*(X(:,2:n)*diag(((1:n-1)./n).^(thet-1))-X(:,1:n-1)*diag(((0:n-2)./n).^(thet-1)))*(b(n-1:-1:1)-b(n:-1:2))';

    f3 = -1/2*(A'*A)*X(:,n);
    f4 = A'*y;

    f = f1+f3+f4;
    x = B\f;
    X(:,n+1) = x;

    Err_199(n+1) = norm(x-x_real,2)/norm(x_real,2);
    Residual(n) = norm(A*x - y);
    x_199 = x;

    if n>10 && Residual(n)>Residual(n-1)
        fprintf('dont decrease')
        % break;
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
    title({'debias \theta = 1.99'});
    figure(3); plot(max(1,n-100):n+1,Err_199(max(1,n-100):n+1))
    title({'L^2 error'});

end
