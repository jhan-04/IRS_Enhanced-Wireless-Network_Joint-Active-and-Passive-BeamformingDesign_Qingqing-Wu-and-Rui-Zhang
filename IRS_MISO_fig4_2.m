%centralized algorithm 함수로 만들어서 fig4 그린것
%+distributed 
%마지막
clear all
clc
i=0;
N=50;
for d=15:5:50
    
    i=i+1;
    %SNR_up(i)=SNR_Upperbound(N,d);
    SNR_cen(i)=SNR_opt_IRS(N,d);%centralized 22
    SNR_dis(i)=SNR_IRS_dis(N,d);%distributed33   
    SNR_AuMRT(i)=SNR_A_u_MRT(N,d);%AP-user MRT444
    SNR_AIMRT(i)=SNR_A_I_MRT(N,d);%5555
    SNR_noIRS(i)=SNR_opt_noIRS(N,d);%%%6
    
end

hold on

%plot([15:5:50],SNR_up,'-o')%1
plot([15:5:50],SNR_cen,'g')%2
plot([15:5:50],SNR_dis,'--b')%3
plot([15:5:50],SNR_AuMRT,'-.^r')%4
plot([15:5:50],SNR_AIMRT,'-.vc')%5
plot([15:5:50],SNR_noIRS,'s:k')%6

grid on
xlabel('AP-user horizontal distance, d');
ylabel('Receive SNR (dB)');
title('Fig. 4 in the Paper');
legend('joint design centralized','joint design distributed','AP-user MRT','AP-IRS MRT','Wirhout IRS');
hold off


