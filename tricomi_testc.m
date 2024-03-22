% Best version so far (d, e, f - have problems)

% tricomi u_xx + xu_yy = f, u = g on boundary

% starting from Neuberger 461 Fall '15: NONhomog Lapalace's Eq. w/NONhomog BC

% x, y, and mode discritizations for plots, BC; Fourier basis matrix
c       = -pi/2; 
d       =  pi/2;

% a       = -2; 
% bump function billiards 

a = -1.69; a=-2.87; a=-3.6; a = -4.55; a = -5.82; a=-7.1; a=-8.25; a = -9.32;
% a = -9.3405;

    
% a = -3.5058;a=-3.51 % w75
%     a=-4.14;
%     a=-6.77;
%     a = -4.28284;
%     
%     a = -1.97; %a = -2.2; %a = -2.4; %a = -4.4; 
%     a=-4.250; a = -4.34; a=-4.32;
    
b       =  1;
n       = 10; 
dy      = (d-c)/n;
m       = floor((b-a)/dy);
N       = n-1; M = m-1;       


dx      = dy;     % change to dx = (b-a)/m;
b       = a+m*(dx); % delete this line

x       = linspace(a,a+m*(dx),m+1)';  % replace a+... with b
xs      = x(2:end-1);         
y       = linspace(c,d,n+1)'; 
ys      = y(2:end-1);         

[Y,X]   = meshgrid(y,x);      
[YS,XS] = meshgrid(ys,xs);    

%         g2=0
%     --------->y
%     |        |
% g1  |        | g3
%     |________|
%     V   g4=0
%     x
%

g1 = @(x)   zeros(size(x));
g3 = @(x)   zeros(size(x)); 
ex = 0;

% g1 = @(x) cos(3*x);
% g3 = @(x) cos(3*x)*airy(3^(2/3)*pi/2)/airy(-3^(2/3)*pi/2);
% exact   = @(x,y) airy(3^(2/3)*x) .* cos(3*y) / airy(-3^(2/3)*pi/2);

% g1 = @(x) cos(x);
% ab      = [airy(a) airy(2,a); airy(b) airy(2, b)] \ [ 1;0];
% exact   = @(x,y) (ab(1)*airy(x)+ab(2)*airy(2,x)).*cos(y); 

% g1 = @(x) cos(3*x);
% ab      = [airy(3^(2/3)*a) airy(2,3^(2/3)*a); airy(3^(2/3)*b) airy(2,3^(2/3)*b)] \ [ 1;0];

% Bump propogated: use a = -1.75; n = 200; b>=0;
s = 3;
g1 = @(x) 1*cos(s*x).*(x<pi/(2*s)).*(x>-pi/(2*s));

% ab      = [airy(-s^(2/3)*pi/2) airy(2,-s^(2/3)*pi/2); airy(s^(2/3)*pi/2) airy(2,s^(2/3)*pi/2)] \ [ 1;0];
% exact   = @(x,y) -(ab(1)*airy(s^(2/3)*x)+ab(2)*airy(2,s^(2/3)*x)).*cos(s*y); 

% g1 = @(x) -4/pi^2*(x-pi/2).*(x+pi/2);

ex = 1;
g1 = @(x) cos(x);
g3 = @(x) cos(3*x);
ab1      = [airy(3^(2/3)*a) airy(2,3^(2/3)*a); airy(3^(2/3)*b) airy(2,3^(2/3)*b)] \ [ 0;1];
ab2      = [airy(a) airy(2,a); airy(b) airy(2,b)] \ [ 1;0];
exact    = @(x,y) (ab1(1)*airy(3^(2/3)*x) + ab1(2)*airy(2, 3^(2/3)*x)) .* cos(3*y) + ...
                  (ab2(1)*airy(x)         + ab2(2)*airy(2, x))         .* cos(  y); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tricomi matrix for square 
In      = speye(N);
R1      = spdiags(xs, 0, M,M);
T1       = spdiags(kron([1 -2 1],ones(M,1)),-1:1,M,M); 
T2       = spdiags(kron([1 -2 1],ones(N,1)),-1:1,N,N); 
L       = kron(In,T1)/dx^2 + kron(T2,R1)/dy^2;

% Initialize RHS to zeros subtract boundary terms g_i defined below   

q        = 3; 
[Ky,Yk]  =  meshgrid(2*(1:q)-1,ys);
V        = cos(Ky.*Yk);
ag       = 2/pi*V'*g1(ys)*dy;
ug       = V * ag;

if ex == 0
    fprintf('using efnct coeff sys solves\n');
    sum = 0;
    for i =1:q
        s = 2*i-1;
        ab  = pinv([airy(s^(2/3)*a) airy(2,s^(2/3)*a); airy(s^(2/3)*b) airy(2,s^(2/3)*b)]) * [ ag(i);0];
        sum = sum + (ab(1)*airy(s^(2/3)*XS)+ab(2)*airy(2,s^(2/3)*XS)).*cos(s*YS); 
    end
else
    fprintf('using provided exact sol\n');
    sum = exact(XS,YS);
end

R        = zeros(M,N);
R(1,:)   = R(1,:)   - ug' / dy^2;          
% R(1,:)   = R(1,:)   - g1(ys') / dy^2;          
R(end,:) = R(end,:) - g3(ys')/dy^2;         
R        = reshape(R,N*M,1);

% Solve system
u        = L \ R;

% add BC for plot and L2-norm comparison with exact solution
u        = [g1(y)'; [zeros(M,1), full(reshape(u, M,N)), zeros(M,1)]; g3(y)'];
v        = [g1(y)'; [zeros(M,1), sum,                   zeros(M,1)]; g3(y)'];

figure(10);
% figure(floor(abs(a))); 
% subplot(2,1,1); 
% contourf(X, Y, u,'EdgeColor','none'); 
subplot(2,1,1); surf(X, Y, u,'EdgeColor','none'); 
subplot(2,1,2); surf(X, Y, v-u,'EdgeColor','none'); 

l2err = norm(reshape(u-v,(n+1)*(m+1),1))/sqrt(n+1)/sqrt(m+1); % R-sum for (\int(u-v)^2 dx dy)^.5 
fprintf('||Exact Sol - D2 Sol|| = %f\n', l2err);

% figure(3); plot(x,u(:,n/2),x,v(:,n/2));

% figure(4);
% W = @(a,b,k) airy(k^(2/3)*a) * airy(2,k^(2/3)*b) - airy(2,k^(2/3)*a) * airy(k^(2/3)*b);
% x5 = linspace(-5,0,500)';
% plot(x5,W(x5,1,1),x5,W(x5,1,3),'.');hold on;
% plot(-1.96729,0,'ko',-2.40983,0,'ko',x5,zeros(size(x5)),'k'); hold off;
