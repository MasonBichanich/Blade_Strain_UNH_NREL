function [cross_corr,lags] = f_cross_corr(x,y,nlags,dt)
%For positive lags, x is leading y

for jj=0:nlags
    sx=0;
    sy=0;
    sxx=0;
    syy=0;
    sxy=0;
    Ng=0;
    for n=1:length(x)-jj
        x1=x(n);
        y1=y(n+jj);
        if isnan(x1) == false && isnan(y1) == false
            Ng=Ng+1;
            sx=sx+x1;
            sy=sy+y1;
            sxx=sxx+x1^2;
            syy=syy+y1^2;
            sxy=sxy+x1*y1;
        end
    end
    cross_corr(nlags+1+jj)=(sxy/Ng-(sx/Ng)*(sy/Ng))/((sxx/Ng-sx^2/Ng^2)^0.5*(syy/Ng-sy^2/Ng^2)^0.5);
end

for jj=1:nlags
    sx=0;
    sy=0;
    sxx=0;
    syy=0;
    sxy=0;
    Ng=0;
    for n=1:length(x)-jj
        x1=x(n+jj);
        y1=y(n);
        if isnan(x1) == false && isnan(y1) == false
            Ng=Ng+1;
            sx=sx+x1;
            sy=sy+y1;
            sxx=sxx+x1^2;
            syy=syy+y1^2;
            sxy=sxy+x1*y1;
        end
    end
    cross_corr(nlags+1-jj)=(sxy/Ng-(sx/Ng)*(sy/Ng))/((sxx/Ng-sx^2/Ng^2)^0.5*(syy/Ng-sy^2/Ng^2)^0.5);
end

lags=-nlags*dt:dt:nlags*dt;

end

