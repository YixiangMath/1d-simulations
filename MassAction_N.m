


I=1;
for a=0.2:0.01:1.2

m = 0;
%x = [0 0.005 0.01 0.05 0.1 0.2 0.5 0.7 0.9 0.95 0.99 0.995 1];
x=0:0.02:1;
T=800;
%t = [0 0.005 0.01 0.05 0.1 0.5 1 1.5 2];
t=0:1:T;

pdex4ic=@(x)pdex4ica(a,x);
sol = pdepe(m,@pdex4pde,pdex4ic,@pdex4bc,x,t);
u1 = sol(:,:,1);
u2 = sol(:,:,2);


%

% plot(x, u1(31,:),'LineWidth',2)
% hold on
% 
% plot(x, u2(31,:),'LineWidth',2)

deltax=1/(length(x)-1);
Total(I)=(sum(u2(T+1,:))-u2(T+1,1)/2-u2(T+1, length(x))/2)*deltax;

I=I+1;
end

h=figure
plot(2.2:0.01:3.2,Total, LineWidth=2) 
xlabel('N');
% ylabel('\int_0^1 I(x, 40)dx');

%xlim([3.5 4.5])
%ylim([0 0.5])
% 
 legend('\int_0^1 I(x, 40)dx')
% 
 saveas(h, '/Users/ywu/Library/Mobile Documents/com~apple~CloudDocs/My Papers/2023-3-degenerate SIS/1d simulations/massN', 'epsc');
% 
% 
% hold off






% --------------------------------------------------------------
function [c,f,s] = pdex4pde(x,t,u,DuDx)
c = [1; 1];                                  
f = [0; 1] .* DuDx;                   
%F=0.5*(1+x).*u(1).*u(2)-(1+x).*(4-pi*sin(pi*x)).*u(2);
F=0.5*(1+x).*u(1).*u(2)-(4-pi*sin(pi*x)).*u(2);
%y = u(1) - u(2);
%F = exp(5.73*y)-exp(-11.47*y);
s = [-F; F];
end 

% --------------------------------------------------------------
function u0 = pdex4ica(a,x);
u0 = [2+cos(pi*x); a+cos(pi*x)];   
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


