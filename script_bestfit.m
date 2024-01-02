clear all, close all, clc

%% Plotting
fig_path = './Figures/';

if ~exist('Figures', 'dir')
    mkdir('Figures')
end

%% Dimensional parameters
R = 1e-4; % radius
N = 150; % number of terms in series expansion
P = 5e-8; % mass transfer coefficient
r = linspace(0,R,101)'; %501 % plot concentration at these positions
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

%% Best-fit diffusivity
options = optimset('Display','iter');
alpha = 80; % heterogeneous D(r)
[mah,tmh] = FGM_model(R,P,Dmin,Dmax,Davg,kmin,kmax,kavg,c0min,c0max,c0avg,r,tc,ts,tm,alpha,N,AbsTol);

alpha = 1e-4; % homogeneous D
F = @(D) sum((mah - FGM_model(R,P,Dmin,Dmax,D,kmin,kmax,kavg,c0min,c0max,c0avg,r,tc,ts,tm,alpha,N,AbsTol)).^2);
D = 5e-13; % initial guess
D = fminsearch(F,D,options)
[mahh] = FGM_model(R,P,Dmin,Dmax,D,kmin,kmax,kavg,c0min,c0max,c0avg,r,tc,ts,tm,alpha,N,AbsTol);

alpha = 1e4; % two-layer D
F = @(Dv) sum((mah - FGM_model(R,P,Dv(1),Dv(2),Davg,kmin,kmax,kavg,c0min,c0max,c0avg,r,tc,ts,tm,alpha,N,AbsTol)).^2);
Dv = [5e-13,5e-12]; % initial guess
Dv = fminsearch(F,Dv,options)
[maht] = FGM_model(R,P,Dv(1),Dv(2),Davg,kmin,kmax,kavg,c0min,c0max,c0avg,r,tc,ts,tm,alpha,N,AbsTol);

%% Plots
figure;
differenceh = mah-mahh;
differencet = mah-maht;
plot(tmh,differenceh,'-','Color','b','LineWidth',line_width)
hold on
plot(tmh,differencet,'-','Color','r','LineWidth',line_width)
set(gca,'Fontsize',font_size,'FontName','Times','Color',background_color)
xlabel('$\hat{t}$','Interpreter','LaTeX')
ylabel('Difference','Interpreter','LaTeX')
ylim([-0.03,0.01])
xlim([0,tmh(end)])
text(14,-0.01,'Homogeneous','Fontsize',font_size,'Color','b','FontName','Times');
text(14,-0.0145,'Two-layer','Fontsize',font_size,'Color','r','FontName','Times');

exportgraphics(gcf,[fig_path,'difference_bestfit.pdf'])