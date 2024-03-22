n = 100;
x = linspace(0,1,n+1)';
[Y,X] = meshgrid(x,x);

f = @(x,y) sin(pi*x) .* sin(pi*y);

Z = f(X,Y);

subplot(2,1,1);contourf(X,Y,Z,'Edgecolor','none');
subplot(2,1,2);contourf(X,Y,Z);