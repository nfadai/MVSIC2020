clear
close all;

%% INPUT

% NS = 65, 129, 257, 513
NR = 1; 
NS = 129;

% q = 1;
sc = 5;
sav = 1;

for q = 1
D = ['q = ',num2str(q)];
disp(D)



if q == 1
Theta=[3 3 8 3 3 .1 .1 8 .5 .1 1 .05 .1 .2 0.05];
elseif q == 2
Theta=[4 4 6 3 3 .1 .1 8 .5 .1 1 .05 .1 .2 0.05];    
elseif q == 3
Theta=[5 5 8 3 3 .1 .1 8 .5 .1 1 .05 .1 .2 0.05];  
end

gam1=Theta(1); %3,4,5
gam2=Theta(2); %3,4,5
kap= Theta(3); %8,6,8
gam3=Theta(4); %3
gam4=Theta(5); %3;
lam1=Theta(6); %.1;
lam2=Theta(7); %.1;
zeta=Theta(8); %8;
rho=Theta(9);%.5;
sig=Theta(10);%.1;
mu=Theta(11);%1;
a=Theta(12);%.05;
b=Theta(13);%.1;
w=Theta(14);%.2;
V0=Theta(15); % 0.05
Y0 = [1,0,V0,0,0];

% efast_var={'$\gamma_1$','$\gamma_2$','$\kappa$','$\gamma_3$','$\gamma_4$','$\lambda_1$','$\lambda_2$','$\zeta$'...
%             '$\rho$','$\sigma$','$\mu$','$\alpha$','$\beta$','$\omega$','$V_0$','dummy'};
        
efast_var={'$\alpha$','$\beta$','$\gamma_1$','$\gamma_2$','$\gamma_3$','$\gamma_4$','$\zeta$','$\kappa$','$\lambda_1$','$\lambda_2$',...
           '$\mu$','$\rho$','$\sigma$','$\omega$','$V_0$','dummy'};
        
        
        
        
pmin=(1/sc)*[Theta,1]; 
pmax=sc*[Theta,1]; 
S = length(pmin);


%% Output

wantedN=NS*S*NR;    % wanted no. of sample points
MI = 4; %: maximum number of fourier coefficients that may be retained in calculating the partial variances without interferences between the assigned frequencies

% OUTPUT
% SI[] : first order sensitivity indices
% STI[] : total effect sensitivity indices
% Other used variables/constants:
% OM[] : vector of k frequencies
% OMi : frequency for the group of interest
% OMCI[] : set of freq. used for the compl. group
% X[] : parameter combination rank matrix
% AC[],BC[]: fourier coefficients
% FI[] : random phase shift
% V : total output variance (for each curve)
% VI : partial var. of par. i (for each curve)
% VCI : part. var. of the compl. set of par...
% AV : total variance in the time domain
% AVI : partial variance of par. i
% AVCI : part. var. of the compl. set of par.
% Y[] : model output

%% PARAMETERS AND ODE SETTINGS (they are included in the following file)

tspan = linspace(0,50,1e3);
time_points= 0:1:50;

% Computation of the frequency for the group of interest OMi and the # of sample points NS (here N=NS)
OMi = floor(((wantedN/NR)-1)/(2*MI)/S);
NS = 2*MI*OMi+1;
if(NS*NR < 65)
    fprintf(['Error: sample size must be >= ' ...
    '65 per factor.\n']);
    return;
end

%% Pre-allocation of the output matrix Y, Y will save only the points of interest specified in the vector time_points
Ya(NS,length(time_points),length(Y0),length(pmin),NR)=0;  % pre-allocation

tic
for i=1:S
    disp(i)
    OMci = SETFREQ(S,OMi/2/MI,i);   
    for L=1:NR
        cj = 1;
        for j=1:S
            if(j==i)
                OM(i) = OMi;
            else
                OM(j) = OMci(cj);
                cj = cj+1;
            end
        end

        FI = rand(1,S)*2*pi; 
        S_VEC = pi*(2*(1:NS)-NS-1)/NS;
        OM_VEC = OM(1:S);
        FI_MAT = FI(ones(NS,1),1:S)';
        ANGLE = OM_VEC'*S_VEC+FI_MAT;
        
        X(:,:,i,L) = 0.5+asin(sin(ANGLE'))/pi; 
        X(:,:,i,L) = parameterdist(X(:,:,i,L),pmax,pmin,0,1,NS,'unif'); 
        
        for run_num=1:NS
            [i run_num L];           
%               [t,y]=ode15s(@(t,y)ODE_efast(t,y,X(:,:,i,L),Theta,run_num),tspan,Y0,[]);   
                
%              [t,y] = ode_solv(X(:,:,i,L),run_num,tspan,Theta,a);            
%             Y1(run_num,:,:,i,L)= output_func(y,t,1);
            
            Y1(run_num,:,:,i,L)= output_func(X(:,:,i,L),run_num,tspan,1);
            Y2(run_num,:,:,i,L)= output_func(X(:,:,i,L),run_num,tspan,2);
            Y3(run_num,:,:,i,L)= output_func(X(:,:,i,L),run_num,tspan,3);
            Y4(run_num,:,:,i,L)= output_func(X(:,:,i,L),run_num,tspan,4);          
            
        end 
    end 
end 
toc

% save Model_efast.mat;
% [Si,Sti,rangeSi,rangeSti] = efast_sd(Y1,OMi,MI,1,1);
% mybarplot(Si,Sti,efast_var)

[Si1,Sti1,rangeSi1,rangeSti1] = efast_sd(Y1,OMi,MI,1,1);
[Si2,Sti2,rangeSi2,rangeSti2] = efast_sd(Y2,OMi,MI,1,1);
[Si3,Sti3,rangeSi3,rangeSti3] = efast_sd(Y3,OMi,MI,1,1);
[Si4,Sti4,rangeSi4,rangeSti4] = efast_sd(Y4,OMi,MI,1,1);

[Si1a, Sti1a] = reorder(Si1, Sti1);
[Si2a, Sti2a] = reorder(Si2, Sti2);
[Si3a, Sti3a] = reorder(Si3, Sti3);
[Si4a, Sti4a] = reorder(Si4, Sti4);

Si = [Si1a,Si2a,Si3a,Si4a];
Sti = [Sti1a,Sti2a,Sti3a,Sti4a];


%%
close all
mybarplot(Si,Sti,efast_var)

if sav == 1
filename = ['New_q',num2str(q),'_sc',num2str(sc),'_NR',num2str(NR),'_NS',num2str(NS)];
save(filename)
end

end


mybarplot(Si,Sti,efast_var)
