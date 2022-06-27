%centralized algorithm
%fig 4



clc
clear all
P=10^(5/10)/1000;%5dBm
P_noise=10^(-80/10)/1000;%-80dBm
M=8;N=50;
d0=51;
dv=2;
D=0;

% 
% Hr_H=(randn(1,N)+1i*randn(1,N))/sqrt(2);%Rayleigh
% Hd_H=(randn(1,M)+1i*randn(1,M))/sqrt(2);%Rayleigh

for i=1:N
    Hr_H(1,i)=1;%/sqrt(2)+1i/sqrt(2);    
end

for i=1:M
    Hd_H(1,i)=1;%/sqrt(2)+1i/sqrt(2);
end


kdb=1000;
k=10^(kdb/10);
g=sqrt(k/(k+1))+sqrt(1/(k+1)).*(randn(1,M)+1i*randn(1,M))/sqrt(2);%Rician

%g=zeros(1,M)+1;

for j=1:N
%20*log10(51) + 20*log10(28*10^6)-147.55=35.5446dB
%g2=g.*sqrt(10^(-3).*10^(5/10).*10^(-35.5446/10));   
%g2=g.*sqrt(10^(-3).*10^(5/10)./(4*pi/10^(5/10)*d0^2)); 
g2=g.*sqrt(10^(-3).*10^(5/10).*d0^-2.2); 
G(j,:)=g2;
end

for d=15:5:50

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



% cvx_begin
% variable V(N+1,N+1) symmetric
% maximize(trace(R*V));
% subject to
% diag(V)==1;
% V == semidefinite(N+1);
% cvx_end
%%%%%%%%%%%%


[U,W,Z] = svds(V);%V=U*W*Z'

obj_max=-100;
 n=length(W);
for j=1:1000
    
    r=sqrt(1/2)*(randn(n,1)+1i*randn(n,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
    v_=U*sqrt(W)*r;
    obj=(v_'*R*v_);
    
    if obj>=obj_max
        obj_max=obj;
        v_max=v_;
    end
    
    
end



t=v_max(N+1);
for j=1:N
    
v(j,:)=exp(1i*angle(v_max(j)/t));
%v(j,:)=v_max(j)/t;
end

% for i=n+1:N
%     
% v(i,:)=0;
% 
% end



w=sqrt(P).*(v'*Phi+hd_H)'./norm(v'*Phi+hd_H);

obj_final=P*(norm(v'*Phi+hd_H))^2;
obj_final2=(abs((v'*Phi+hd_H)*w))^2;



SNR=10*log10(obj_final/P_noise);%dB
 D=D+1;
SNR_re(D)=SNR;

    
end
hold on
plot([15:5:50],SNR_re)%d~SNR



grid on
xlabel('AP-user horizontal distance, d');
ylabel('Receive SNR (dB)');
title('Fig. 4 in the Paper');
axis([])
hold off








