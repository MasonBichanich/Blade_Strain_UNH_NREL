function [ S,f ] = f_spectra( ts,dt,nbands,nens )

N=length(ts);

%ensemble average

K=floor(N/nens);
jj=0:K/2;
fk=jj/K/dt;
df=(fk>0 & fk<1/2/dt);
fk=fk(df);

So=0;
for k=1:nens
    tsk=ts((k-1)*K+1:k*K);
    tsk=tsk-mean(tsk);
    Gj=(1/K)*fft(tsk);
    G_pro=Gj.*conj(Gj);
    Sj(:,k)=2*K*dt*G_pro;
    So=So+Sj(df,k);

end
    
Sk=So/nens;

%band average
B=floor(length(Sk)/nbands);

for b=1:B
    S(b)=mean(Sk((b-1)*nbands+1:b*nbands));
    f(b)=mean(fk((b-1)*nbands+1:b*nbands));
end

end 


