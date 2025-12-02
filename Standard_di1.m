%%%%%%%This program solves u(x, t) and graph the solution
%%%%Divide x=[0, 1] and t=[0, T] into even intervals%%%%
clc
clear
dx=0.02;
ds=1; di=0;
dt=0.0002;
%%%%Solve the equation in the time interval [0, T]
T=40;
%T=8;

%%%M, N is the number of nodes corresponding to space and time 
M=1/dx+1;
N=T/dt+1;

%%%%%%%%%U stores the solution u(x, t)
%S=zeros(N,M);
%I=zeros(N,M);
%%%%%%Initial Condition%%%%%%%%%%%%%%%%%%%
for J=2:(M-1)
    %U(1,I)=sin(pi*(I-1)*dx);
    %S(1,I)=2*(I-1)*dx*(1- (I-1)*dx);
    S(1,J)=2+cos(pi*(J-1)*dx);
    I(1,J)=1.5+cos(pi*(J-1)*dx);
end

%%%%%%%%Boundary on the initial time t=0%%%%%%%%%%%%
S(1,1)=S(1,2);
S(1,M)=S(1,M-1);
I(1,1)=I(1,2);
I(1,M)=I(1,M-1);


%%%%%%%%%Finite difference scheme%%%%%%%%%
for K=2:N     %%%%%%%%%Time loop
    for J=2:(M-1) %%Space loop

        S(K,J)=S(K-1, J)+ dt/(dx)^2*ds*( S(K-1,J-1)-2*S(K-1, J) +S(K-1, J+1))-dt*(2- (abs((J-1)*dx-0.5))^0.5)*S(K-1,J)*I(K-1,J)/(S(K-1,J)+I(K-1,J)) +dt*1.5* I(K-1,J); 

        I(K,J)=I(K-1, J)+ dt/(dx)^2*di*( I(K-1,J-1)-2*I(K-1, J) +I(K-1, J+1))+dt*(2-(abs((J-1)*dx-0.5))^0.5)*S(K-1,J)*I(K-1,J)/(S(K-1,J)+I(K-1,J) ) -dt*1.5* I(K-1,J);     

        %%%%%%%%Calculate boundary values 
        S(K,1)=S(K,2);
        S(K,M)=S(K,M-1);
        I(K,1)=I(K,2);
        I(K,M)=I(K,M-1);
    end
end

%%%%%%%%Graph the the solution at T
h=figure
plot(0:dx:1, S(N,:),'LineWidth',2)
hold on
plot(0:dx:1, I(N,:),'LineWidth',2)

xlabel('x')
%ylabel('Population density')
legend('S(x, 40)', 'I(x, 40)')

 saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-3-degenerate SIS/1d simulations/standardDI1', 'epsc');
