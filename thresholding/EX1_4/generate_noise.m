
m = 100; tau = 1/m; s = (0:tau:1)'; t = 0:tau:1;
K = min(s,t).*(1-max(s,t));
t1 = (0:tau:1/5)'; t2 = (1/5+tau:tau:2/5)';
t3 = (2/5+tau:tau:3/5)';
t4 = (3/5+tau:tau:4/5)';t5 = (4/5+tau:tau:1)';
x_real = [zeros(size(t1));-ones(size(t2));zeros(size(t3));ones(size(t4));zeros(size(t5));];
A = tau*K*diag([1/2,ones(1,m-1),1/2]);% Rectangle integral formula
y_real = A*x_real;
noise = (2*rand(size(y_real))-1);save('noise','noise');