function out = output_func(X,run_num,tspan,opt)

gam1= X(run_num,1); %3,4,5
gam2=X(run_num,2); %3,4,5
kap= X(run_num,3); %8,6,8
gam3=X(run_num,4); %3
gam4=X(run_num,5); %3;
lam1=X(run_num,6); %.1;
lam2=X(run_num,7); %.1;
zeta=X(run_num,8); %8;
rho=X(run_num,9);%.5;
sig=X(run_num,10);%.1;
mu=X(run_num,11);%1;
a=X(run_num,12);%.05;
b=X(run_num,13);%.1;
w=X(run_num,14);%.2;
V0=X(run_num,15); % 0.05
dummy=X(run_num,16);

n=2;
m=2;

zeta*kap/(gam2+kap);

e=@(t) 0.*(t>40);
phi=@(t) 0*(1-exp(20-t)).*(t>20);

f=@(v) (v).^n./(b^n+(v).^n);
g=@(v) (v).^m./(w^m+(v).^m);

S=@(t,x) [x(1).*(1-x(1)-kap*(1-e(t))*x(3));
    kap*(1-e(t)).*x(1).*x(3)-gam1.*x(2)*(1+lam1*x(4));
    zeta*x(2)-gam2*x(3).*(1+lam2*x(4))-kap*(1-e(t))*x(1)*x(3);
    gam3*(x(5)-x(4))+rho*x(1).*x(3)./(a+x(3));
    (x(1)+sig*x(2)).*f(x(3))+mu*x(4).*g(x(3))-gam4*(1+phi(t)).*x(5)];


options=odeset('RelTol',3e-14,'absTol',1e-16,'NonNegative',1:5);

sol = ode15s(@(t,x) S(t,x),tspan,[1 0 V0 0.0 0],options); %0.05

Y=sol.y;

if opt == 1
    %first output metric: min of S
    out = min(Y(1,:));

elseif opt == 2 
    %average amount of fluctuating cytokine (>0)
    T=max(tspan);
%    dCdt =  (Y(1,:)+sig*Y(2,:)).*f(Y(3,:))+...
%        mu*Y(4,:).*g(Y(3,:))-gam4*(1+phi(sol.x)).*Y(5,:);
%    ds = sqrt(1+(dCdt/max(Y(5,:))).^2);
    out =1/(max(Y(5,:))*T)*trapz(sol.x(2:end),abs(diff(Y(5,:))./diff(sol.x)));
    
elseif opt == 3
    %max-min cytokine difference after transient response (t=20)    
    [~,Ti]=min(abs(sol.x-20));
    out = (max(Y(5,Ti:end))-min(Y(5,Ti:end)))/max(Y(5,:));
    
elseif opt == 4
    %total accrued cytokines after transient response (t=20)   
    [~,Ti]=min(abs(sol.x-20));
    out = trapz(sol.x(Ti:end),Y(5,Ti:end))/max(Y(5,:));
end

end