function [V, lam, s] = SL_pqn(p,q,u,x,k)

% solves (py')' + qy = -\lambda uy, 0-D BC

% INPUT: 
% p,q,u = functions defining SL problem, assuming 0-D,
% x = n+1 length pt grid for [a,b], 
% k = # eigs to compute

% OUTPUT:
% V = cols of evcts, 
% lam = evals, 
% s = <e,e> = \int e^2 u dx = +/-1 (weighted semi inner prod)

n = length(x)-1;            % # divisions = length cell grid 
a = x(1); b = x(end); dx = x(2)-x(1);
v = x(1:end-1) + dx/2;      % cell grid
X = kron(v, ones(1,k));

Dinn = spdiags(kron([-1 1],ones(n+1,1)),-1:0, n+1, n)   / dx;
Dinn(1,1)=0; Dinn(end,end)=0;
Dout = spdiags(kron([-1 1],ones(n+1,1)), 0:1, n,   n+1) / dx;
P    = spdiags(p(x),    0, n+1, n+1);
Q    = spdiags(q(v),    0, n,   n);
UI   = spdiags(1./u(v), 0, n,   n);
L    = UI * (Dout * P * Dinn + Q);

% norm eigs w/weighted semi-inner prod, order by eval mag, cannonical flip
[V,D]       = eigs(L+speye(size(L)),k,'sm');
lam         = diag(D)'-1;
[~, ndx]    = sort(abs(lam));
lam         = lam(ndx);
V           = V(:,ndx);
snrm        = diag(V'*(u(X).*V)*dx);
s           = sign(snrm);
nrm         = sqrt(abs(snrm));
V           = V * diag(sign(V(1,:))'./nrm); 

end
