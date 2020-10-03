function mybarplot(Si,Sti)

efast_var_hor={'$\alpha$','$\beta$','$\gamma_1$','$\gamma_2$','$\gamma_3$','$\gamma_4$','$\zeta$','$\kappa$','$\lambda_1$','$\lambda_2$',...
           '$\mu$','$\rho$','$\sigma$','$\omega$','$V_0$','dum'};
        

l = size(Si,1);
x = 1:l;
    
for j = 1:size(Si,2)

figure
y1 = Si(:,j);
y2 = Sti(:,j);
y = [y1, y2];
h = bar(x, y);
set(gcf,'Position', [10 10 800 500]);
legend(h,{'First Order','Total Order'})
ylim([0 1])
xticks(1:16)
xticklabels(efast_var_hor)
ylabel('eFAST Sensitivity','Fontsize',22')
xlabel('Parameter','Fontsize',22')
shg

end


end