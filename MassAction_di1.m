




m = 0;
%x = [0 0.005 0.01 0.05 0.1 0.2 0.5 0.7 0.9 0.95 0.99 0.995 1];
x=0:0.02:1;
T=50;
%t = [0 0.005 0.01 0.05 0.1 0.5 1 1.5 2];
t=0:1:T;

sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t);
u1 = sol(:,:,1);
u2 = sol(:,:,2);


h=figure

plot(x, u1(51,:),'LineWidth',2)
hold on

plot(x, u2(51,:),'LineWidth',2)

xlabel('x')

legend('S(x,50)', 'I(x,50)')

saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-3-degenerate SIS/1d simulations/massDI1', 'epsc');


hold off




% h=figure
% 
% plot(x, u1(101,:),'LineWidth',2)
% hold on
% 
% plot(x, u2(101,:),'LineWidth',2)
% 
% xlabel('x')
% 
% legend('S(x,100)', 'I(x,100)')
% 
% saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-3-degenerate SIS/1d simulations/stDS1T1000', 'epsc');
% 
% 
% hold off
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
f = [1; 0] .* DuDx;                   
F=u(1).*u(2)-(4-pi*sin(pi*x)).*u(2);
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


% %%%%%%%This program solves u(x, t) and graph the solution
% %%%%Divide x=[0, 1] and t=[0, T] into even intervals%%%%
% clc
% clear
% dx=0.02;
% ds=1; di=0;
% dt=0.0002;
% %%%%Solve the equation in the time interval [0, T]
% T=10;
% %T=8;
% 
% %%%M, N is the number of nodes corresponding to space and time 
% M=1/dx+1;
% N=T/dt+1;
% 
% %%%%%%%%%U stores the solution u(x, t)
% %S=zeros(N,M);
% %I=zeros(N,M);
% %%%%%%Initial Condition%%%%%%%%%%%%%%%%%%%
% for J=2:(M-1)
%     %U(1,I)=sin(pi*(I-1)*dx);
%     %S(1,I)=2*(I-1)*dx*(1- (I-1)*dx);
%     S(1,J)=2+cos(pi*(J-1)*dx);
%     I(1,J)=1.5+cos(pi*(J-1)*dx);
% end
% 
% %%%%%%%%Boundary on the initial time t=0%%%%%%%%%%%%
% S(1,1)=S(1,2);
% S(1,M)=S(1,M-1);
% I(1,1)=I(1,2);
% I(1,M)=I(1,M-1);
% 
% 
% %%%%%%%%%Finite difference scheme%%%%%%%%%
% for K=2:N     %%%%%%%%%Time loop
%     for J=2:(M-1) %%Space loop
% 
%         S(K,J)=S(K-1, J)+ dt/(dx)^2*ds*( S(K-1,J-1)-2*S(K-1, J) +S(K-1, J+1))-dt*1*S(K-1,J)*I(K-1,J) +dt*( 4-pi*sin(pi* (J-1)*dx) )* I(K-1,J); 
% 
%         I(K,J)=I(K-1, J)+ dt/(dx)^2*di*( I(K-1,J-1)-2*I(K-1, J) +I(K-1, J+1))+dt*1*S(K-1,J)*I(K-1,J) -dt*( 4-pi*sin(pi* (J-1)*dx) )* I(K-1,J);     
% 
%         %%%%%%%%Calculate boundary values 
%         S(K,1)=S(K,2);
%         S(K,M)=S(K,M-1);
%         I(K,1)=I(K,2);
%         I(K,M)=I(K,M-1);
%     end
% end
% 
% %%%%%%%%Graph the the solution at T
% h=figure
% plot(0:dx:1, S(N,:),'LineWidth',2)
% hold on
% plot(0:dx:1, I(N,:),'LineWidth',2)
% 
% xlabel('x')
% %ylabel('Population density')
% legend('S', 'I')
% 
%  saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-3-degenerate SIS/1d simulations/massDI1', 'epsc');
% % subplot(1,2,1)
% %  [X,Y]=meshgrid( 0:dx:1,0:dt:T);
% % mesh(X,Y,S)
% % xlabel('x')
% % ylabel('t')
% % zlabel('S')
% % subplot(1,2,2)
% %  [X,Y]=meshgrid( 0:dx:1,0:dt:T);
% % mesh(X,Y,I)
% % xlabel('x')
% % ylabel('t')
% % zlabel('I')
% % sum=0;
% % for K=2:N
% %     sum(K)=0;
% %     for J=1:M
% %     sum(K)=sum(K)+dx*S(K,M)+dx*I(K,M);
% %     end
% % end