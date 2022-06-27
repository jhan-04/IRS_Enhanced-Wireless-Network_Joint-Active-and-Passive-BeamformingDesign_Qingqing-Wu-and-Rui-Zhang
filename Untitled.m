%centralized algorithm function222
N=50;
d=43;
num=0;

for c=1:1000
N
d
c
    P=10^(5/10)/1000;%5dBm
    P_noise=10^(-80/10)/1000;%-80dBm
    M=8;
    d0=51;
    dv=2;
    
    Hr_H=(randn(1,N)+1i*randn(1,N))/sqrt(2);%Rayleigh
    Hd_H=(randn(1,M)+1i*randn(1,M))/sqrt(2);%Rayleigh
    
    
    
    kdb=100;
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
    
    
    Phi=diag(hr_H)*G;
    R_=[Phi*Phi' Phi*hd;hd_H*Phi' 0];
    %     R=[Phi*Phi' Phi*hd;hd_H*Phi' 0];
    
    r=round(log10(mean(mean(R_))));
    R=R_;
   % R=round(R_,real(6-r));
    ishermitian(R);








[U,W,Z] = svd(R);%U==Z

v_max=W(1,1)*U(:,1);







    
   
    
    t=v_max(N+1);
    for j=1:N
        
        v(j,:)=exp(1i*angle(v_max(j)/t));%angle of reflection
        %v(j,:)=v_max(j)/t;
    end
    theta=diag(v');
    

    w2=sqrt(P).*(hr_H*theta*G+hd_H)'./norm(hr_H*theta*G+hd_H);
    w=sqrt(P).*(v'*Phi+hd_H)'./norm(v'*Phi+hd_H);
    
    obj_final=P*(norm(v'*Phi+hd_H))^2;
    obj_final2=(abs((v'*Phi+hd_H)*w))^2;
    value_new=(abs((hr_H*theta*G+hd_H)*w2))^2;
    
    SNR_c(c)=10*log10(obj_final/P_noise);%dB
    obj_final3(c)= obj_final;
    %   end
    
    
end
SNR=mean(SNR_c);
%SNR=10*log(mean(obj_final3)/P_noise);







