clear all
close all
dbstop if error
%parameter in the jump distribution:
k=15;
%M here is nsim...
nsim=10000000;

iw=0;%if unity, work with Weibull, otherwise, mixture of normals

%parameters of the Weibull distribution
a=10;
b=20;
%parameters of mixture of Normals distribution
sig1=.02;%standard deviation of first normal
sig2=.01;%standard deviation of second normal
mu1=-.06;%mean of first
mu2=.06;%mean of second
mix=.5;%probability of normal 1

if iw == 1
    xx=7:.001:11.5;
    Y = wblpdf(xx,a,b);
    [maxY,i]=max(Y);
    mode=xx(i);
else
    x1=min([mu1 mu2])-5*max([sig1 sig2]);
    x2=max([mu1 mu2])+5*max([sig1 sig2]);
    ngrid=5000;
    inc=(x2-x1)/ngrid;
    xx=x1:inc:x2;
    [Y]=mixednormpdf(xx,mu1,mu2,sig1,sig2,mix);
    [maxY,i]=max(Y);
    mode=xx(i);  
end
    

%Get the second derivative about the mode of the Weibull distribution, 
%for the purpose of obtaining the jump distribution for the MCMC algorithm.
%The second derivative can be approximated numerically as the first derivative
%of the first derivative. A numerical approximation of the function f around
%the point, x, is (f(x+eps)-f(x-eps))/(2*eps), so the approximation of
%the second derivative is (eps1=eps):
%fpp=(f(x+eps+eps1)-f(x-eps+eps1)-(f(x+eps-eps1)-f(x-eps-eps1)))/(2*eps*2*eps1);

eps=.000001;
x=mode;

%get the variance required for the Laplace approximation and the MCMC:
if iw == 1
    lf1=log(wblpdf(x+2*eps,a,b));
    lf2=log(wblpdf(x-2*eps,a,b));
elseif iw == 0
    [Ya]=mixednormpdf(x+2*eps,mu1,mu2,sig1,sig2,mix);
    [Yb]=mixednormpdf(x-2*eps,mu1,mu2,sig1,sig2,mix);
    lf1=log(Ya);
    lf2=log(Yb);
end
lfpp=(lf1-2*log(maxY)+lf2)/(2*eps*2*eps);
sVN=sqrt(-1/lfpp);
sV=sVN;

randn('seed',0);
thet(1)=mode;
ix=0;
tic
for ii = 2:nsim
    rr=randn*sV;
    u=rand;
    x=thet(ii-1)+k*rr;
    x0=thet(ii-1);
    if x <= 0 && iw == 1
        thet(ii)=thet(ii-1);%the Weibull random variable is non-negative, so a negative candidate has density zero.
    else
        if iw == 1
            fn=wblpdf(x,a,b);
            fd=wblpdf(x0,a,b);
        elseif iw == 0
            fn=mixednormpdf(x, mu1,mu2,sig1,sig2,mix);
            fd=mixednormpdf(x0,mu1,mu2,sig1,sig2,mix);
        end
        lam=fn/fd;
        if u < lam
            thet(ii)=x;
            ix=ix+1;
        else
            thet(ii)=thet(ii-1);
        end
    end
end
toc
ia=min([1000,round(.01*nsim)]);
[N,X]=hist(thet([ia:nsim]),150);
b1=max(diff(X));
b2=min(diff(X));
if abs(b1-b2) > .1e-10
    error('(gammaMCMC) something went wrong')
end
bb=b1;

%Next, scale N so that the Riemann integral over the histogram is unity.
%Here is the logic....
%The Riemann integral is the sum of the area of rectangles under the
%integral. In the case of a histogram, it is the sum of the width of
%the rectangles times their heigth. The hist command provides a height, N(i),
%associated with each of a set of 150 bins of width b, and these have to be
%scaled by a constant independent of i, to ensure that the Riemann integral
%is unity. The scaled height is n(i). The selected constant is a.
%(It can be verified numerically that a similar adjustment before graphing
%a pdf function is redundant, i.e., a is approximately unity.)

aa=bb*sum(N);
n=N/aa;
I=find(X<max(xx)&X>min(xx));
s2=max(I);
s1=min(I);

aa=['MCMC, k = ',num2str(k),', % acceptance = ',num2str(100*ix/nsim)];
II=[s1:s2];
plot(xx,Y,X(II),normpdf(X(II),mode,sVN),'x-',X(II),n(II))
if iw == 1
    legend('Weibull','Laplace approximation',aa)
    ff=['Weibull distribution, mode = ',num2str(mode),', a = ',num2str(a),', b = ',num2str(b), ...
        ', number of MCMC simulations = ',num2str(nsim-ia)];
    title(ff,'FontSize',24)
elseif iw == 0
    legend('Mixture of Normals','Laplace approximation',aa)
    ff1=['mixture of normals, mu1 = ',num2str(mu1),', mu2 = ',num2str(mu2),', sig1 = ',num2str(sig1), ...
        ', sig2 = ',num2str(sig2),', mix = ',num2str(mix)];
    ff2=['number of MCMC simulations = ',num2str(nsim-ia)];
    title({ff1;ff2},'FontSize',24)
end
axis tight

1+1;
