function [ S12,f,Coh2,Ph ] = my_cross_spectra( ts1,ts2,dt,nbands,nens )

N=length(ts1);

%ensemble average

K=floor(N/nens);
jj=0:K/2;
fk=jj/K/dt;
df=(fk>0 & fk<1/2/dt);
fk=fk(df);

So1=0;
So2=0;
Co=0;
Qo=0;
for k=1:nens
    tsk=ts1((k-1)*K+1:k*K);
    tsk=tsk-mean(tsk,'omitnan');
    G1=(1/K)*fft(tsk);
    A1=real(G1);
    B1=imag(G1);
    Sj1(:,k)=K*dt*(A1.^2+B1.^2);
    So1=So1+Sj1(df,k);
    
    tsk=ts2((k-1)*K+1:k*K);
    tsk=tsk-mean(tsk,'omitnan');
    G2=(1/K)*fft(tsk);
    A2=real(G2);
    B2=imag(G2);
    Sj2(:,k)=K*dt*(A2.^2+B2.^2);
    So2=So2+Sj2(df,k);
    
    Cj(:,k)=K*dt*(A1.*A2+B1.*B2);
    Co=Cj(df,k)+Co;
    Qj(:,k)=K*dt*(B1.*A2-A1.*B2);
    Qo=Qj(df,k)+Qo;

end
    
Sk1=So1/nens;
Sk2=So2/nens;
Ck=Co/nens;
Qk=Qo/nens;

%band average
B=floor(length(Sk1)/nbands);

for b=1:B
    S1(b)=mean(Sk1((b-1)*nbands+1:b*nbands));
    S2(b)=mean(Sk2((b-1)*nbands+1:b*nbands));
    C(b)=mean(Ck((b-1)*nbands+1:b*nbands));
    Q(b)=mean(Qk((b-1)*nbands+1:b*nbands));
    f(b)=mean(fk((b-1)*nbands+1:b*nbands));
end

S12=C-1i*Q;
Coh2=(C.^2+Q.^2)./S1./S2;
Ph=atan2(-1*Q,C)*180/pi;

end 


