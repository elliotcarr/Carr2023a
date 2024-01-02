clear all, close all, clc

%% Plotting
fig_path = './Figures/';

if ~exist('Figures', 'dir')
    mkdir('Figures')
end

%% Dimensional parameters
R = 1e-4; % radius
P = 5e-8; % mass transfer coefficient
r = linspace(0,R,501)'; %501 % plot concentration at these positions
tc = [10^1,10^2,10^3,10^4]; % plot concentration at these times (1D plot)
ts = [0.05,0.1,0.5]*10^4; % plot concentration at these times (2D plot)
tm = linspace(0,3*10^4,101); % plot mass at these times

% D(r) = Dmax + (Dmin-Dmax)[0.5+atan(alpha*(r-sigma)/R)/pi]
% k(r) = kmin + (kmax-kmin)[0.5+atan(alpha*(r-sigma)/R)/pi]
alpha_vec = [1e-4,20,80,1e4]; % alpha values in D(r), k(r) etc
Dmin = 1e-13; % absolute min diffusivity (as alpha -> infty)
Dmax = 1e-11; % absolute max diffusivity (as alpha -> infty)
kmin = 0; % absolute min reaction rate (as alpha -> infty)
kmax = 0; % absolute max reaction rate (as alpha -> infty)
c0min = 0.4;
c0max = 0.4;
c0avg = 0.4;
AbsTol = 1e-9; %integral tolerance

% plotting options
font_size = 30;
line_width = 3;
background_color = [1,1,1];

% Average values of D and k
Davg = 3/R^3*(Dmax*((R/2)^3)/3 + Dmin*(R^3-(R/2)^3)/3);
kavg = 3/R^3*(kmin*((R/2)^3)/3 + kmax*(R^3-(R/2)^3)/3);

%% Eigenvalue analysis
alpha = 80; % heterogeneous D(r)

% benchmark number of terms in series expansion
Np = 10;
[mah,tmh,chNold] = FGM_model(R,P,Dmin,Dmax,Davg,kmin,kmax,kavg,c0min,c0max,c0avg,r,tc,ts,tm,alpha,Np,AbsTol);
diffN = Inf;
while diffN > 1e-6
    Np = Np + 10;
    [mah,tmh,chN] = FGM_model(R,P,Dmin,Dmax,Davg,kmin,kmax,kavg,c0min,c0max,c0avg,r,tc,ts,tm,alpha,Np,AbsTol);
    Np
    diffN = max(abs(chN - chNold),[],'all');
    chNold = chN;
end

N = Np;
Nvec = 10:10:50; % eigenvalues to test

[mah,tmh,chb,rh] = FGM_model(R,P,Dmin,Dmax,Davg,kmin,kmax,kavg,c0min,c0max,c0avg,r,tc,ts,tm,alpha,N,AbsTol);
err = zeros(size(Nvec));
for i = 1:length(Nvec)
    [mah,tmh,chN] = FGM_model(R,P,Dmin,Dmax,Davg,kmin,kmax,kavg,c0min,c0max,c0avg,r,tc,ts,tm,alpha,Nvec(i),AbsTol);
    err(i) = max(abs(chN-chb),[],'all');
end

%% Plots
figure
semilogy(Nvec,err,'.-','Color','b','LineWidth',2,'MarkerSize',30)
set(gca,'Fontsize',font_size,'FontName','Times','Color',background_color,'YTick',10.^(-[5:-1:0]))
xlabel('$N$','Interpreter','LaTeX')
ylabel('Error','Interpreter','LaTeX')
ylim([0.99e-5,1e0])
xlim([Nvec(1),Nvec(end)])
box on
drawnow

fig1 = figure(1);
exportgraphics(fig1,[fig_path,'eigenvalue_analysis.pdf'])