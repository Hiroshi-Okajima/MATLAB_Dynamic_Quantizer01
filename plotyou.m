function plotyou(te,Ap,Bp,Cp,A,B,C,d)

[u,uq,yq,y,e] = sim(te,Ap,Bp,Cp,A,B,C,d);
max_abs_e = max(abs(e));
max_abs_y = max(abs(y));
fprintf('max|e| = %4.3f  max|y| = %4.3f \n\n',max_abs_e,max_abs_y)

% 理想出力yと量子化後の出力yq
figure(1)
stairs(0 : te,y,'k','linewidth',1.5)
hold on
stairs(0 : te,yq,'r','linewidth',3)
hold off
axis([0 te -6 6])
xlabel('$t$', 'interpreter', 'latex','fontsize',18)
ylabel({'$y(t),y_q(t)$'}, 'interpreter', 'latex','fontsize',18)
legend({'$y(t)$', '$y_q(t)$'}, 'interpreter', 'latex','fontsize',18)
grid on

% 量子化前の入力uと量子化後の入力(量子化器からの出力)uq
figure(2)
stairs(0 : te,u,'k','linewidth',3)
hold on
stairs(0 : te,uq,'r','linewidth',1.5)
hold off
axis([0 te -3.2 3.2])
xlabel('$t$', 'interpreter', 'latex','fontsize',18)
ylabel({'$u(t),u_q(t)$'}, 'interpreter', 'latex','fontsize',18)
legend({'$u(t)$', '$u_q(t)$'}, 'interpreter', 'latex','fontsize',18)
grid on

% 出力誤差e=yq-y
figure(3)
stairs(0 : te,e,'k','linewidth',2)
axis([0 te -0.5 0.5])
xlabel('$t$', 'interpreter', 'latex','fontsize',18)
ylabel({'$e(t)$'}, 'interpreter', 'latex','fontsize',18)
legend({'$e(t)$'}, 'interpreter', 'latex','fontsize',18)
grid on

function [u,uq,yq,y,e] = sim(te,Ap,Bp,Cp,A,B,C,d)

np = length(Ap);
nq = length(A);

u       = zeros(te+1,1);
uq      = zeros(te+1,1);
yq      = zeros(te+1,1);
y       = zeros(te+1,1);

for t = 0 : te
    if(t==0)
        xi = zeros(nq,1); % ξ(0)  
        xq = zeros(np,1); % xq(0)         
        x  = zeros(np,1); %  x(0)                     
    end
 
    u(t+1)     = 0.3*sin(0.6*t)+0.4*sin(0.1*t)+0.3*sin(0.3*t); % u(t) 
    ubar       = C*xi+u(t+1);                     % u^-(t) = Cξ(t)+u(t)
    uq(t+1)    = d*round(1/2+ubar/d)-d/2;    % uq(t)  = Qst[u^-(t)]
    yq(t+1)    = Cp*xq;                           % yq(t)  = Cp*xq(t)
    y(t+1)     = Cp*x;                            % y(t)   = Cp*x(t)
    
    xi = A  * xi + B  * (uq(t+1)-u(t+1));  % ξ(t+1)   =  Aξ(t)+ B*(uq(t)-u(t))
    xq = Ap * xq + Bp *  uq(t+1);          % xq(t+1)   = Apξ(t)+Bp*uq(t)
    x  = Ap *  x + Bp *   u(t+1);          % x(t+1)    = Apξ(t)+Bp* u(t)
end

e = yq-y;
