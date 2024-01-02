function [mah,tmh,ch,rh,mahinf,Dh,k,tsh,csh] = FGM_model(R,P,Dmin,Dmax,Davg,kmin,kmax,kavg,c0min,c0max,c0avg,r,tc,ts,tm,alpha,N,AbsTol)

% Diffusivity
fdf = @(sigma,r) 0.5+atan(alpha*(r-sigma)/R)/pi;
Ddf = @(sigma,r) Dmax + (Dmin-Dmax)*fdf(sigma,r);
% g = @(sigma) integral(@(r)r.^2.*Ddf(sigma,r),0,R) - R^3/3*Davg;
g = @(sigma) integral(@(r)r.^2.*Ddf(sigma,r),0,R)/(R^3/3*Davg) - 1;
%options = optimset('Display','iter');
sigma = fzero(g,R/2);
fprintf(['alpha = ',num2str(alpha),'\n']);
fprintf(['sigma = ',num2str(sigma,'%1.3e'),'\n']);
absg = abs(g(sigma));
fd = @(r) fdf(sigma,r);
D = @(r) Ddf(sigma,r);
fprintf(['avg(D(r)) = ',num2str(3/R^3*integral(@(r) r.^2.*D(r),0,R)),'\n']);
fprintf(['Deff = ',num2str(R^3/(3*integral(@(r) ((r.^2)./D(r)),0,R))),'\n']);

% Reaction rate
k = @(r) kmin + (kmax-kmin)*fd(r);
fprintf(['avg(k(r)) = ',num2str(3/R^3*integral(@(r) r.^2.*k(r),0,R)),'\n']);

% Initial concentration
c0 = @(r) c0max + (c0min-c0max)*fd(r);
fprintf(['avg(c0(r)) = ',num2str(3/R^3*integral(@(r) r.^2.*c0(r),0,R)),'\n\n']);

% Non-dimensional variables (h for hat)
Dh = @(rh) D(rh*R)/Dmax;
kh = @(rh) k(rh*R)*R^2/Dmax;
c0h = @(rh) c0(rh*R)/c0avg;
[Rh,Ph,rh,tch,tsh,tmh,sigmah] = deal(R/R,P*R/Dmax,r/R,Dmax*tc/R^2,Dmax*ts/R^2,Dmax*tm/R^2,sigma/R);

% Eigenvalues and Eigenfunctions
lambda = eigenvalues(Ph/Dh(1),1,Rh,N);
X = @(n,r) Xfunc(n,r,lambda);
Xd = @(n,r) Xdfunc(n,r,lambda);

% Matrix A
syms rs
Dhr = diff(Dh(rs),rs);
Dhr = matlabFunction(Dhr,'Vars',rs);

A = zeros(N,N); T0 = zeros(N,1); B = zeros(N,N);
for n = 1:N
    for m = 1:N
        A(m,n) = integral(@(r) Dhr(r).*Xd(n,r).*X(m,r),0,Rh,'AbsTol',AbsTol) ...
            - integral(@(r) r.^2.*(lambda(n)^2*Dh(r)+kh(r)).*X(n,r).*X(m,r),0,Rh,'AbsTol',AbsTol);
    end
    T0(n) = integral(@(r) r.^2.*c0h(r).*X(n,r),0,Rh);
end

% Concentration (1D plot)
ch = zeros(length(rh),length(tch));
for j = 1:length(tch)
    T = expm(tch(j)*A)*T0;
    for n = 1:N
        ch(:,j) = ch(:,j) + T(n)*X(n,rh);
    end
end

% Concentration (2D plot)
csh = zeros(length(rh),length(tsh));
for j = 1:length(tsh)
    T = expm(tsh(j)*A)*T0;
    for n = 1:N
        csh(:,j) = csh(:,j) + T(n)*X(n,rh);
    end
end

% Released Mass
chR = [0, zeros(1,length(tmh)-1)];
scal = integral(@(r) r.^2.*c0h(r),0,Rh);
for j = 2:length(tmh)
    Ttd = (expm(tmh(j)*A)-eye(size(A)))*(A\T0);
    for n = 1:N
        chR(j) = chR(j) + Ttd(n)*X(n,Rh);
    end
end
mah = Ph*chR/scal;

% Total Released Mass
tinf = 1e6;
Ttd = (expm(tinf*A)-eye(size(A)))*(A\T0);
chRinf = 0;
for n = 1:N
    chRinf = chRinf + Ttd(n)*X(n,Rh);
end
mahinf = Ph*chRinf/scal;

end

function Xr = Xfunc(n,r,lambda)

if r == 0
    Xr = 2*sqrt(lambda(n))/(sqrt(2*lambda(n)-sin(2*lambda(n))))*lambda(n);
else
    Xr = 2*sqrt(lambda(n))/(sqrt(2*lambda(n)-sin(2*lambda(n))))*(sin(lambda(n)*r)./r);
end

end

function Xdr = Xdfunc(n,r,lambda)

Xdr = 2*sqrt(lambda(n))/(sqrt(2*lambda(n)-sin(2*lambda(n))))*(lambda(n)*r.*cos(lambda(n)*r) - sin(lambda(n)*r));

end