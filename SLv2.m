function SLv1

% solves (py')' + qy = -\lambda uy, 0-D BC
% Tested for Airy's SOV

n = 200; k = n-2; a = -1; b = 1; dx = (b-a)/n;
x = linspace(a,b,n+1)'; v = x(1:end-1) + dx/2;
X = kron([a;v;b], ones(1,k)); 

mm = 3;
p = @(x)  ones(size(x));
q = @(x)  zeros(size(x));
u = @(x)  x.^mm; 

Dinn = spdiags(kron([-1 1],ones(n+1,1)),-1:0, n+1, n)   / dx;
Dout = spdiags(kron([-1 1],ones(n+1,1)), 0:1, n,   n+1) / dx;
P    = spdiags(p(x),    0, n+1, n+1);
Q    = spdiags(q(v),    0, n,   n);
UI   = spdiags(1./u(v), 0, n,   n);
L    = UI * (Dout * P * Dinn + Q);

[V,D]       = eigs(-L,k,'sm');
lam         = diag(D)';
[~, ndx]    = sort(abs(lam));
lam         = lam(ndx);
V           = V(:,ndx);
snrm        = diag(V'*(X(2:end-1,:).^mm.*V)*dx);
s           = sign(snrm);
nrm         = sqrt(abs(snrm));
V           = V * diag(1./nrm); 

Y = [zeros(1,k); V ;zeros(1,k)];
figure(1); plot(X(:,1:10),Y(:,1:10)); ylim([-10,10]);

[lam(1:6)]

% Test: expand f by efncts w/weighted signed-norm, compare
f = @(x)  (1/4 - x.^2).*(x>-1/2).*(x<1/2);
% f = @(x)  (x>-1/2).*(x<1/2);


c = s .* (f(v)'*(X(2:end-1,:).^mm.*V*dx))';
w = V*c;

figure(2); plot(v,f(v),'b', v,w,'r');

% % near singular frequencies
% figure(3);
% pp=find(lam<0);qq=find(lam>0);
% nnlam = -lam(pp);
% plot(1./sin(sqrt(nnlam))); 
% 
% figure(11);
% plot(X(:,pp),Y(:,pp));
% 
% figure(12);
% plot(X(:,qq),Y(:,qq));
end

