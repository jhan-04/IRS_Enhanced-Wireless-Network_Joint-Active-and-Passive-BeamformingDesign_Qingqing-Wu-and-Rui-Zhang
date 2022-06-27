%centralized algorithm함수 만들어서 fig5
%d=15,43,50일때 N(the number of refelcting element)가 증가하면 어떻게 SNR이 바뀌는지


clear all
clc


figure(1)
 d=15;
    i=0;
    for N=30:10:100
        i=i+1;
        %SNR_up(i)=SNR_Upperbound(N,d);
        SNR_cen(i)=SNR_opt_IRS(N,d);%centralized 22
        SNR_dis(i)=SNR_IRS_dis(N,d);%distributed33
        SNR_AuMRT(i)=SNR_A_u_MRT(N,d);%AP-user MRT444
        SNR_AIMRT(i)=SNR_A_I_MRT(N,d);%5555
        SNR_noIRS(i)=SNR_opt_noIRS(N,d);%%%6
        
        
    end
    
    
    hold on
    
    %plot([30:10:100],SNR_up,'-o')%1
    plot([30:10:100],SNR_cen,'g')%2
    plot([30:10:100], SNR_dis,'--b')%3
    plot([30:10:100],SNR_AuMRT,'-.^r')%4
    plot([30:10:100],SNR_AIMRT,'-.vc')%5
    plot([30:10:100],SNR_noIRS,'s:k')%6


grid on
xlabel('Number of reflecting elements,N');
ylabel('Receive SNR (dB)');
title('Fig. 5 in the Paper d=15m');

legend('joint design centralized','joint design distributed','AP-user MRT','AP-IRS MRT','Wirhout IRS');
%legend('d=15 with IRS','d=15 without IRS','d=50 with IRS','d=50 without IRS');
hold off


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
 d=50;
    i=0;
    for N=30:10:100
        i=i+1;
        %SNR_up(i)=SNR_Upperbound(N,d);
        SNR_cen(i)=SNR_opt_IRS(N,d);%centralized 22
        SNR_dis(i)=SNR_IRS_dis(N,d);%distributed33
        SNR_AuMRT(i)=SNR_A_u_MRT(N,d);%AP-user MRT444
        SNR_AIMRT(i)=SNR_A_I_MRT(N,d);%5555
        SNR_noIRS(i)=SNR_opt_noIRS(N,d);%%%6
        
        
    end
    
    
    hold on
    
    %plot([30:10:100],SNR_up,'-o')%1
    plot([30:10:100],SNR_cen,'g')%2
    plot([30:10:100], SNR_dis,'--b')%3
    plot([30:10:100],SNR_AuMRT,'-.^r')%4
    plot([30:10:100],SNR_AIMRT,'-.vc')%5
    plot([30:10:100],SNR_noIRS,'s:k')%6


grid on
xlabel('Number of reflecting elements,N');
ylabel('Receive SNR (dB)');
title('Fig. 5 in the Paper d=50');

legend('joint design centralized','joint design distributed','AP-user MRT','AP-IRS MRT','Wirhout IRS');
%legend('d=15 with IRS','d=15 without IRS','d=50 with IRS','d=50 without IRS');
hold off








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
figure(3)
 
d=43;
i=0;
for N=30:10:100
    
    i=i+1;
    SNR_cen(i)=SNR_opt_IRS(N,d);%centralized 22
    SNR_dis(i)=SNR_IRS_dis(N,d);%distributed33
    SNR_AuMRT(i)=SNR_A_u_MRT(N,d);%AP-user MRT444
    SNR_AIMRT(i)=SNR_A_I_MRT(N,d);%5555
    SNR_noIRS(i)=SNR_opt_noIRS(N,d);%%%6
    
end
hold on
plot([30:10:100],SNR_cen,'g')%2
plot([30:10:100], SNR_dis,'--b')%3
plot([30:10:100],SNR_AuMRT,'-.^r')%4
plot([30:10:100],SNR_AIMRT,'-.vc')%5
plot([30:10:100],SNR_noIRS,'s:k')%6
grid on
xlabel('Number of reflecting elements,N');
ylabel('Receive SNR (dB)');
title('Fig. 5 in the Paper d=43');
legend('joint design centralized','joint design distributed','AP-user MRT','AP-IRS MRT','Wirhout IRS');


