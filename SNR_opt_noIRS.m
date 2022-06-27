%IRS없을경우 SNR666
function SNR=SNR_opt_noIRS(N,d)

P=10^(5/10)/1000;%5dBm
P_noise=10^(-80/10)/1000;%-80dBm
M=8;
d0=51;
dv=2;
for c=1:2000
 

Hd_H=(randn(1,M)+1i*randn(1,M))*sqrt(0.5);%Rayleigh

for i=1:M
    Hd_H(1,i)=1/sqrt(2)+1i/sqrt(2);
    Hd_H(1,i)=sqrt(0.5)*(-1)^round(rand()*10,0)+1i*sqrt(0.5)*(-1)^round(rand()*10,0);  
end



d_hd=sqrt(d^2+dv^2);%d1

hd_H=Hd_H.*sqrt(10^(-3).*10^(-10/10).*d_hd^-3);
hd=hd_H';




obj_final= P*(norm(hd_H))^2;



SNR_c(c)=10*log10(obj_final/P_noise);


end

SNR=mean(SNR_c);


end










