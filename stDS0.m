

m = 0;
%x = [0 0.005 0.01 0.05 0.1 0.2 0.5 0.7 0.9 0.95 0.99 0.995 1];
x=0:0.02:1;
T=100;
%t = [0 0.005 0.01 0.05 0.1 0.5 1 1.5 2];
t=0:1:T;

sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t);
u1 = sol(:,:,1);
u2 = sol(:,:,2);


h=figure

plot(x, u1(101,:),'LineWidth',2)
hold on

plot(x, u2(101,:),'LineWidth',2)

xlabel('x')

legend('S(x,100)', 'I(x,100)')

saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-3-degenerate SIS/1d simulations/stDS0T100', 'epsc');


hold off



% 
% h=figure
% 
% plot(x, u1(1001,:),'LineWidth',2)
% hold on
% 
% plot(x, u2(1001,:),'LineWidth',2)
% 
% xlabel('x')
% 
% legend('S(x,1000)', 'I(x,1000)')
% 
% saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-3-degenerate SIS/1d simulations/stDS1T1000', 'epsc');
% 
% 
% hold off
% 
% 
% 
% 
% 
% 
% 
% h=figure
% 
% plot(x, u1(T,:),'LineWidth',2)
% hold on
% 
% plot(x, u2(T,:),'LineWidth',2)
% 
% xlabel('x')
% 
% legend('S(x,20000)', 'I(x,20000)')
% 
% saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-3-degenerate SIS/1d simulations/stDS1Large', 'epsc');
% 
% 
% hold off

% figure
% surf(x,t,u1)
% title('u1(x,t)')
% xlabel('Distance x')
% ylabel('Time t')
% 
% figure
% surf(x,t,u2)
% title('u2(x,t)')
% xlabel('Distance x')
% ylabel('Time t')
% --------------------------------------------------------------
function [c,f,s] = pdex4pde(x,t,u,DuDx)
c = [1; 1];                                  
f = [0; 1] .* DuDx;                   
F=(1+sin(pi*x)).*u(1).*u(2)./(u(1)+u(2))-1.5*u(2);
%y = u(1) - u(2);
%F = exp(5.73*y)-exp(-11.47*y);
s = [-F; F];
end 

% --------------------------------------------------------------
function u0 = pdex4ic(x);
u0 = [2+cos(pi*x); 1.5+cos(pi*x)];   
end

% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t)
%pl = [0; ul(2)];  
pl = [0; 0]; 
%ql = [1; 0];
ql = [1; 1];
%pr = [ur(1)-1; 0];     
pr = [0; 0]; 
qr = [1; 1]; 
end
