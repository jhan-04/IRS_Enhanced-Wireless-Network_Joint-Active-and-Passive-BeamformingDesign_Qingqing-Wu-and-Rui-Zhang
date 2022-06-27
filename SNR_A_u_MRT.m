%funtion of distributed algorithm444

function SNR=SNR_A_u_MRT(N,d)


for c=1:1000
    
P=10^(5/10)/1000;%5dBm
P_noise=10^(-80/10)/1000;%-80dBm
M=8;
d0=51;
dv=2;

Hr_H=(randn(1,N)+1i*randn(1,N))/sqrt(2);%Rayleigh
Hd_H=(randn(1,M)+1i*randn(1,M))/sqrt(2);%Rayleigh
% for i=1:N
%     Hr_H(1,i)=sqrt(0.5)*(-1)^round(rand()*10,0)+1i*sqrt(0.5)*(-1)^round(rand()*10,0);  
% end
% 
% for i=1:M
%     Hd_H(1,i)=sqrt(0.5)*(-1)^round(rand()*10,0)+1i*sqrt(0.5)*(-1)^round(rand()*10,0);
% end




% for i=1:N
%     Hr_H(1,i)=sqrt(0.5)+1i*sqrt(0.5);
% end
% 
% for i=1:M
%     Hd_H(1,i)=sqrt(0.5)+1i*sqrt(0.5);
% end
% 
% 







kdb=1000;
k=10^(kdb/10);
g=sqrt(k/(k+1))+sqrt(1/(k+1)).*(randn(1,M)+1i*randn(1,M))/sqrt(2);%Rician
for j=1:N
g2=g.*sqrt(10^(-3).*10^(5/10).*d0^-2.2); 
G(j,:)=g2;
end


d_hr=sqrt((d0-d)^2+dv^2);%d2
d_hd=sqrt(d^2+dv^2);%d1

hr_H=Hr_H.*sqrt(10^(-3).*10^(-10/10).*10^(5/10).*d_hr^-3);
hr=hr_H';
hd_H=Hd_H.*sqrt(10^(-3).*10^(-10/10).*d_hd^-3);
hd=hd_H';



w_k=sqrt(P).*hd/norm(hd);


Phase0=angle(hd_H*w_k);

for n=1:1:N
    th_k(n)=Phase0-angle(hr_H(n))-angle(G(n,:)*w_k);%1XN
end
v_H=exp(1i.*th_k);
v=v_H';%Theta=diag(v_H)

% w_MRT=sqrt(P).*hd/norm(hd);
% a=-angle(hd_H*w_MRT);
% w_k=w_MRT.*exp(a);





value_new=(abs((hr_H*diag(v_H)*G+hd_H)*w_k))^2;


SNR_c(c)=10*log10(value_new/P_noise);
end

SNR=mean(SNR_c);

end

