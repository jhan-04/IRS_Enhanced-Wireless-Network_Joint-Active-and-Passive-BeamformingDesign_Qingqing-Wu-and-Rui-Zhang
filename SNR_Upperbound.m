%centralized algorithm function1111결과이상함..


function SNR=SNR_Upperbound(N,d)

P=10^(5/10)/1000;%5dBm
P_noise=10^(-80/10)/1000;%-80dBm
M=8;
d0=51;
dv=2;

% Hr_H=(randn(1,N)+1i*randn(1,N))/sqrt(2);%Rayleigh
% Hd_H=(randn(1,M)+1i*randn(1,M))/sqrt(2);%Rayleigh

for i=1:N
    Hr_H(1,i)=1/sqrt(2)+1i/sqrt(2);    
end

for i=1:M
    Hd_H(1,i)=1/sqrt(2)+1i/sqrt(2);
end


kdb=1000;
k=10^(kdb/10);
g=sqrt(k/(k+1))+sqrt(1/(k+1)).*(randn(1,M)+1i*randn(1,M))/sqrt(2);%Rician

for j=1:N
%20*log10(51) + 20*log10(28*10^6)-147.55=35.5446dB
%g2=g.*sqrt(10^(-3).*10^(5/10).*10^(-35.5446/10));   
g2=g.*sqrt(10^(-3).*10^(5/10).*d0^-2.2); 
G(j,:)=g2;
end


d_hr=sqrt((d0-d)^2+dv^2);%d2
d_hd=sqrt(d^2+dv^2);%d1


hr_H=Hr_H.*sqrt(10^(-3).*10^(-10/10).*10^(5/10).*d_hr^-3);
hr=hr_H';
hd_H=Hd_H.*sqrt(10^(-3).*10^(-10/10).*d_hd^-3);
hd=hd_H';
% hr_H=Hr_H.*10^(-1.5).*10^(-5/10).*d_hr^-1.5.*10^(2.5/10);
% hr=hr_H';
% hd_H=Hd_H.*10^(-1.5).*10^(-5/10).*d_hd^-1.5;
% hd=hd_H';



Phi=diag(hr_H)*G;
R_=[Phi*Phi' Phi*hd;hd_H*Phi' 0];
R=round(R_,15);
ishermitian(R)

%%%%%%%%%%%%
cvx_begin
variable V(N+1,N+1) symmetric
maximize(trace(R*V))
subject to
for n=1:N+1
   V(n,n)==1;
end
V == semidefinite(N+1);
cvx_end

%%%%%%%%%%%%
obj_final=P*trace(R*V);

SNR=10*log10(obj_final/P_noise);%dB
 
end










