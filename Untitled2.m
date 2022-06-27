%centralized algorithm function222
N=50;
d=40;
num=0;
N
for c=1:100
    
    P=10^(5/10)/1000;%5dBm
    P_noise=10^(-80/10)/1000;%-80dBm
    M=8;
    d0=51;
    dv=2;
    
    Hr_H=(randn(1,N)+1i*randn(1,N))/sqrt(2);%Rayleigh
    Hd_H=(randn(1,M)+1i*randn(1,M))/sqrt(2);%Rayleigh
    
    
    
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
    
    
    Phi=diag(hr_H)*G;
    R_=[Phi*Phi' Phi*hd;hd_H*Phi' 0];
    %     R=[Phi*Phi' Phi*hd;hd_H*Phi' 0];
    
    r=round(log10(mean(mean(R_))));
    R=R_;
   % R=round(R_,real(6-r));
    ishermitian(R);
    
%     %%%%%%%%%%%%
%     cvx_begin
%     variable V(N+1,N+1)  hermitian
%     maximize(real(trace(R*V)));
%     subject to
%     for n=1:N+1
%         V(n,n)==1;
%     end
%     V == semidefinite(N+1);
%     cvx_end
%     %%%%%%%%%%%%


%rank(V);


I1=eye(N+1,N+1);
I1(1,1)=1;

I2=eye(N+1,N+1);
I2(2,2)=1;


    %%%%%%%%%%%%
    cvx_begin
    variable yi(N+1,1)
    maximize sum(yi);
    subject to

    R-diag(yi) == semidefinite(N+1);
    cvx_end
    %%%%%%%%%%%%


V=R-diag(yi);

    
%    %%%%%%%%%%%%
%     cvx_begin
%     variable V(N+1,N+1)  hermitian
%     maximize real((R*V)-diag(yi)*V);
%     subject to
%      V == semidefinite(N+1);
%     cvx_end
%     %%%%%%%%%%%%
%     
    
    
    
    
    

    
    
    
    %[U,W,Z] = svds(V);%V=U*W*Z'
    [U,W,Z] = svd(V);
    obj_max=-1000;
    n=length(W);
    v_max=0;
    
    [A,B,C]=svd(R*V);

    Z=0;
    for i=1:N+1
        z=U(:,i)*sqrt(W(i,i));
        Z=Z+z*z';
        
        v_=z;
        amp=abs(v_);
        v_=v_./amp;
        obj=(v_'*R_*v_);
        
        if obj>=obj_max
            obj_max=obj;
            v_max=v_;
        end
        
    end    
    

    
    
    
    
    
    
    
    
    
    
    
  
%       %%%%%%%%%%%%%%%%%%%%%%%
%     
%     %  app = normrnd(zeros(n,1),V);
%     mu=zeros(1+N,1);
%     sigma=V;
%     
%     for i=1:5000
%         
%         sample = mvnrnd(mu,sigma) ;
%         v_=sample.'; 
%         mag=v_'*v_;
%         v_=v_/sqrt(mag);
%        
%         obj=(v_'*R_*v_);
%         if obj>=obj_max
%             obj_max=obj;
%             v_max=v_;
%         end
%         
%         
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%
    
    
    

    
    
    %%%%%%%
   for i=1:n
        r=zeros(n,1);
        r(i,1)=1;
        
        v_=U*W^(1/2)*r;
        amp=abs(v_);
        v_=v_./amp;
        obj=(v_'*R_*v_);
        
        if obj>=obj_max
            obj_max=obj;
            v_max=v_;
        end
        
    end
    %%%%%%%%%
    
    


%     %%%%
% 
%     for j=1:1000
%         
%         r=sqrt(0.5)*(randn(n,1)+1i*randn(n,1));%sqrt(var/2)*(randn(1,N)+1i*randn(1,N))
%         %r=(rand(n,1)+1i*rand(n,1));
%         v_=U*W^(1/2)*r;
%        % v_=V^(1/2)*A*r;
%         amp=abs(v_);
%         v_=v_./amp;
%         obj=(v_'*R*v_);
%         if obj>=obj_max
%             obj_max=obj;
%             v_max=v_;
%         end
%         
%     end
%     
%     %%%%%%%
    
    
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








