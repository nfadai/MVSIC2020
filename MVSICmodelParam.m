function out=MVSICmodelParam(X,opt,PHI)

if nargin==0
    
% This vector of parameter values corresponds to the severe inflammation 
% case (Section 2). Other parameter values in Table 2 are commented next to
% their respective parameter values.

    X=[3 3 8 3 3 .1 .1 8 .5 .1 1 .05 .1 .2 0.05];
    
% eFAST output choice
    opt=3;
end
set(0, 'DefaultLineLineWidth', 2);
gam1= X(1); %3,4,5
gam2=X(2); %3,4,5
kap= X(3); %8,6,8


gam3=X(4); %3
gam4=X(5); %3;
lam1=X(6); %.1;
lam2=X(7); %.1;
zeta=X(8);%8;

rho=X(9);%.5;
sig=X(10);%.1;
mu=X(11);%1;


a=X(12);%.05;
b=X(13);%.1;
w=X(14);%.2;
n=2;
m=2;



gamma1-zeta*kap/(gam2+kap); %if negative, infectious steady-state exists

e=@(t) 0.*(t>20);

phi=@(t) PHI*(1-exp(20-t)).*(t>20); % PHI corresponds to the increase of cytokine clearance due to ICS at t=20.

f=@(v) (v).^n./(b^n+(v).^n);
g=@(v) (v).^m./(w^m+(v).^m);

S=@(t,x) [x(1).*(1-x(1)-kap*(1-e(t))*x(3));
    kap*(1-e(t)).*x(1).*x(3)-gam1.*x(2)*(1+lam1*x(4));
    zeta*x(2)-gam2*x(3).*(1+lam2*x(4))-kap*(1-e(t))*x(1)*x(3);
    gam3*(x(5)-x(4))+rho*x(1).*x(3)./(a+x(3));
    (x(1)+sig*x(2)).*f(x(3))+mu*x(4).*g(x(3))-gam4*(1+phi(t)).*x(5)];

time=linspace(0,50,1e3);

options=odeset('RelTol',1e-13,'absTol',1e-14,'NonNegative',1:5);

sol = ode15s(@(t,x) S(t,x),time,[1 0 X(15) 0.0 0],options); %0.05

Y=sol.y;

if opt==1
    %first output metric: min of S
    out = min(Y(1,:));
    
elseif opt ==2 
    %average amount of fluctuating cytokine (>0)
    T=max(time);

    out =1/(max(Y(5,:))*T)*trapz(sol.x(2:end),abs(diff(Y(5,:))./diff(sol.x)));
    
elseif opt==3
    %max-min cytokine difference after transient response (t=20)
    
    [~,Ti]=min(abs(sol.x-20));
    
    out = (max(Y(5,Ti:end))-min(Y(5,Ti:end)))/max(Y(5,:));
else
    out=NaN;
end


%% Figure plots: uncomment to generate plots of five species (MVSIC)

%  
%  plot(sol.x,Y(1,:),'k',sol.x,Y(2,:),'g',sol.x,Y(3,:),'b',...
%      sol.x,Y(4,:)/max(Y(4,:)),'c',sol.x,Y(5,:)/max(Y(5,:)),'r')
% set(gca,'FontSize',20)
% xlabel('Rescaled Time','FontSize',20)
% ylabel('Rescaled Concentration','FontSize',20)
% 
% 
% axis([0 max(time) 0 1])
% 
% 
% legend({'Susceptible','Infected','Virus','Immune','Cytokine'},'Location','Northeast')
% box on


