close all
clear all
clc

U = [-1 1];
M = 2;
%M = 2;
L = 100;
%[Ap,Bp,Cp] = ssdata(c2d(tf([1 0.5],[1 3 2]),0.05));
%[Ap,Bp,Cp] = ssdata(c2d(tf([1],[1 3 2]),0.1)); % rei1
[Ap,Bp,Cp] = ssdata(c2d(tf([1 20],[1 3 2]),0.1)); % rei2
%[Ap,Bp,Cp] = ssdata(c2d(tf([1 -5],[1 3 2]),0.1)); % rei3

te     = 150;
k_max  = 300;   
m      = 1000;

pa = 4;

xg = zeros(k_max+1,pa);
fg = zeros(k_max+1,1);

x     = zeros(m,pa);
v     = zeros(m,pa);
fxp   = zeros(m,1); 
fx    = zeros(m,1); 
    
a = [-1 1];    
for p = 1 : pa  
    for i = 1 : m  
        x(i,p)  = a(1)+(a(2)-a(1))*rand; 
        v(i,p)  = a(1)+(a(2)-a(1))*rand;                                 
    end
end

xp = x;
    
for k = 0 : k_max
    for i = 1 : m  
        A = [   0      1;
             x(i,1) x(i,2)];
        B = [0;
             1];
        C = [x(i,3) x(i,4)];
        [~,f] = d_f_calculation(L,A,B,C,U(2),U(1),M,Ap,Bp,Cp,te,1,1);         
        fx(i,1)  = f;
        if(k==0)
            fxp(i,1) = fx(i,1);
        else
            if fx(i,1) < fxp(i,1)
                xp(i,:)  = x(i,:); 
                fxp(i,1) = fx(i,1);
            end
        end
    end

    [fg(k+1,1), ig] = min(fxp);
    xg(k+1,:)     = xp(ig,:);
    A = [   0      1;
         xg(k+1,1) xg(k+1,2)];
    B = [0;
         1];
    C = [xg(k+1,3) xg(k+1,4)];
    [d,~] = d_f_calculation(L,A,B,C,U(2),U(1),M,Ap,Bp,Cp,te,1,0);    
    fprintf('k=%d/%d d=%4.4f fg=%4.4f\n',k,k_max,d,fg(k+1,1))
    for i = 1 : m
        for p = 1 : pa
            v(i,p) = 0.9*v(i,p)+1*rand*(xp(i,p)-x(i,p))+1*rand*(xg(k+1,p)-x(i,p)); 
        end
            x(i,:) = x(i,:) + v(i,:);
    end
end

A = [   0      1;
     xg(k_max+1,1) xg(k_max+1,2)];
B = [0;
     1];
 
C = [xg(k_max+1,3) xg(k_max+1,4)];

[d,f] = d_f_calculation(L,A,B,C,U(2),U(1),M,Ap,Bp,Cp,te,1,1);

fprintf('d = %4.4f\nE(Q) = %4.4f\n',d,f)
plotyou(te,Ap,Bp,Cp,A,B,C,d)

figure(100)
plot(0:k_max , fg)
axis([-k_max*0.05 k_max 0.3 3])
xlabel('$k$', 'interpreter', 'latex','fontsize',18)
ylabel({'$E(Q)$'}, 'interpreter', 'latex','fontsize',18)
grid on
hold off

