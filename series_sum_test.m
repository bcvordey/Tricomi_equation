x = linspace(-pi,2*pi,1000)';
M = 100; kk = 0:M;

[X,K]=meshgrid(x,kk);
y = sum(X.^(2*K)./factorial(2*K).*(-1).^K);

x11 = linspace(x(1),x(end),11);
plot(x,y,x11,cos(x11),'r*');