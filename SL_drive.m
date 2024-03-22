
n = 200; k = n-2; a = -1; b = 1; dx = (b-a)/n;
x = linspace(a,b,n+1)'; v = x(1:end-1) + dx/2; xv = [a;v;b];
X = kron(v, ones(1,k)); 

mm = 1;
p = @(x)  ones(size(x));
q = @(x)  zeros(size(x));
u = @(x)  x.^mm; 

[V, lam, s] = SL_pqu(p,q,u,x,k);

j = 1;
e = V(:,j); 
figure(12); plot(xv,[0;e;0]); fprintf('lam(%d) = %f\n',j,lam(j));

% % Test: expand f by efncts w/weighted signed-norm, compare
f = @(x)  (1/4 - x.^2).*(x>-1/2).*(x<1/2);
% f = @(x)  (x>-1/2).*(x<1/2);

c = s .* (f(v)'*(u(X).*V)*dx)';

w = V*c;

figure(2); plot(v,f(v),'b', v,w,'r');
% 
% % % near singular frequencies
% % figure(3);
% % pp=find(lam<0);qq=find(lam>0);
% % nnlam = -lam(pp);
% % plot(1./sin(sqrt(nnlam))); 
% % 
% % figure(11);
% % plot(X(:,pp),Y(:,pp));
% % 
% % figure(12);
% % plot(X(:,qq),Y(:,qq));

