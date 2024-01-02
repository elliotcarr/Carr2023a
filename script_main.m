clear all, close all, clc

%% Plotting
fig_path = './Figures/';

if ~exist('Figures', 'dir')
    mkdir('Figures')
end

% cap_plots = true; % plot 3d capsule plots (true/false)
cap_plots = false;

%% Test Cases
case_name = '1'; % pure diffusion (produces plots in Figs 3, 4, 5, 8, 9 from paper)
% case_name = '2'; % reaction diffusion (produces plots in Figs 3, 6, 7, 10 from paper)

%% Dimensional parameters
R = 1e-4; % radius
N = 100; % number of terms in series expansion
P = 5e-8; % mass transfer coefficient
r = linspace(0,R,501)'; % plot concentration at these positions
tc = [10^1,10^2,10^3,10^4]; % plot concentration at these times (1D plot)
ts = [0.05,0.1,0.5]*10^4; % plot concentration at these times (2D plot)
tm = linspace(0,3*10^4,101); % plot mass at these times

% D(r) = Dmax + (Dmin-Dmax)[0.5+atan(alpha*(r-sigma)/R)/pi]
% k(r) = kmin + (kmax-kmin)[0.5+atan(alpha*(r-sigma)/R)/pi]
alpha_vec = [1e-4,20,80,1e4]; % alpha values in D(r), k(r) etc
Dmin = 1e-13; % absolute min diffusivity (as alpha -> infty)
Dmax = 1e-11; % absolute max diffusivity (as alpha -> infty)
if isequal(case_name,'1')
    kmin = 0; % absolute min reaction rate (as alpha -> infty)
    kmax = 0; % absolute max reaction rate (as alpha -> infty)
    c0min = 0.4;
    c0max = 0.4;
elseif isequal(case_name,'2')
    kmin = 0.8e-4; % absolute min reaction rate (as alpha -> infty)
    kmax = 1e-4; % absolute max reaction rate (as alpha -> infty)
    c0min = 0.4;
    c0max = 0.4;
end
c0avg = 0.4;
AbsTol = 1e-9; %integral tolerance

% Average values of D and k
Davg = 3/R^3*(Dmax*((R/2)^3)/3 + Dmin*(R^3-(R/2)^3)/3);
kavg = 3/R^3*(kmin*((R/2)^3)/3 + kmax*(R^3-(R/2)^3)/3);

% Check for presence of colormaps
Nc = 128; % number of colors
if ~isfile('plasma.m')
    warning(['Download plasma colormap from https://www.mathworks.com/matlabcentral/fileexchange',...
        '/62729-matplotlib-perceptually-uniform-colormaps and place in current directory.'])
    cmap = parula(Nc); % colormap  
else
    cmap = plasma(Nc); % colormap
end

% plotting options
font_size = 30;
line_width = 3;
colors = [0,0,0; 1,0,0; 0.2,0.8,0.2; 0,0,1];
background_color = [1,1,1];
linestyle = {'-','-','-','-'};
label_array = {'(a)','(b)','(c)','(d)'};

for i = 1:length(alpha_vec)

    alpha = alpha_vec(i);
    [mah,tmh,ch,rh,mahinf,Dh,k,tsh,csh] = FGM_model(R,P,Dmin,Dmax,Davg,kmin,kmax,kavg,c0min,c0max,c0avg,r,tc,ts,tm,alpha,N,AbsTol);

    %% Figures
    figure(1); % Concentration
    for j = 1:length(tc)
        leg(1) = plot(rh,ch(:,j),linestyle{i},'Color',colors(i,:),'LineWidth',line_width);
        hold on
    end
    ylim([0,ceil(c0max/c0avg)])
    set(gca,'Fontsize',font_size,'FontName','Times','Color',background_color)
    xlabel('$\hat{r}$','Interpreter','LaTeX')%('r/R')
    ylabel('$\hat{c}(\hat{r},\hat{t})$','Interpreter','LaTeX')%('c(r,t)/C_{0}')
    box on
    text(0.01,0.3,['{ }$\alpha = ',num2str(alpha_vec(i)),'$'],'Units','normalized',...
        'Color',colors(i,:),'Interpreter','LaTeX','FontSize',font_size)
    text(-0.2,-0.17,[label_array{i}],'Units','normalized',...
        'Color','k','FontSize',font_size+12)
    annotation('arrow',[0.85,0.6],[0.85,0.6],'linewidth',line_width,...
        'headlength',15,'headwidth',15)
    drawnow
    pause(1)
    exportgraphics(gcf,[fig_path,'conc',case_name,num2str(i),'.pdf'])
    hold off

    figure(2); % Mass    
    hold on
    plot(tmh,mah,linestyle{i},'Color',colors(i,:),'LineWidth',line_width)
    rt = interp1(mah,tmh,0.99*mahinf); % release time
    plot([rt,rt],[0,0.99*mahinf],'--','Color',colors(i,:),'LineWidth',1)
    set(gca,'Fontsize',font_size,'FontName','Times','Color',background_color)
    xlabel('$\hat{t}$','Interpreter','LaTeX')
    ylabel('$\hat{M}(\hat{t})$','Interpreter','LaTeX')
    ylim([0,1])
    xlim([0,tmh(end)])
    box on
    if i == 4
        text(1.01,0.9,['{ }$\alpha = ',num2str(alpha_vec(1)),'$'],'Units','normalized',...
            'Color',colors(1,:),'Interpreter','LaTeX','FontSize',font_size)
        text(1.01,0.8,['{ }$\alpha = ',num2str(alpha_vec(2)),'$'],'Units','normalized',...
            'Color',colors(2,:),'Interpreter','LaTeX','FontSize',font_size)
        text(1.01,0.7,['{ }$\alpha = ',num2str(alpha_vec(3)),'$'],'Units','normalized',...
            'Color',colors(3,:),'Interpreter','LaTeX','FontSize',font_size)
        text(1.01,0.6,['{ }$\alpha = ',num2str(alpha_vec(4)),'$'],'Units','normalized',...
            'Color',colors(4,:),'Interpreter','LaTeX','FontSize',font_size)
        dar = get(gca,'Position');
        pos = get(gcf,'Position');
        set(gcf,'Position',[pos(1:2),pos(3)*1.24,pos(4)])
        set(gca,'Position',[dar(1:2),dar(3)*1/1.24,dar(4)]);
        drawnow
    end
    if isequal(case_name,'1')
        figure(3); % Diffusivity
        plot(rh,Dh(rh),linestyle{i},'Color',colors(i,:),'LineWidth',line_width)
        hold on
        if i == 1
            plot(rh,Dmin/Dmax*ones(size(rh)),'--','Color',[0.7,0.7,0.7],'LineWidth',line_width)
            plot(rh,Dmax/Dmax*ones(size(rh)),'--','Color',[0.7,0.7,0.7],'LineWidth',line_width)
        end
        ylim([Dmin/Dmax-0.05*(Dmax/Dmax-Dmin/Dmax),1.05*Dmax/Dmax])
        xlabel('$r$','Interpreter','LaTeX')
        set(gca,'Fontsize',font_size,'FontName','Times','Color',background_color,...
            'Xtick',[0,0.5,1],'YTick',[0.01,Davg/Dmax,1],'Layer','top','TickLabelInterpreter','LaTeX')
        labely = ylabel('$D(r)$','Interpreter','LaTeX');
        labely.Position(1) = -0.051; 
        labely.Position(2) = 0.5; 
        if i == 1
            text(-0.17,-0.17,'(a)','Units','normalized','Color','k','FontSize',font_size+12)
        end
        %set(gca,'YTickLabel',{'$D_{\rm{min}}$','$D_{\rm{avg}}$','$D_{\rm{max}}$'}','TickLabelInterpreter','LaTeX')
        set(gca,'YTickLabel',{'$D_{\rm min}$','$D_{\rm avg}$','$D_{\rm max}$'})
        set(gca,'XTickLabel',{'0','$R/2$','$R$'}')
        if i == 4
            text(0.6,0.7,['{ }$\alpha = ',num2str(alpha_vec(1)),'$'],'Units','normalized',...
                'Color',colors(1,:),'Interpreter','LaTeX','FontSize',font_size)
            text(0.6,0.6,['{ }$\alpha = ',num2str(alpha_vec(2)),'$'],'Units','normalized',...
                'Color',colors(2,:),'Interpreter','LaTeX','FontSize',font_size)
            text(0.6,0.5,['{ }$\alpha = ',num2str(alpha_vec(3)),'$'],'Units','normalized',...
                'Color',colors(3,:),'Interpreter','LaTeX','FontSize',font_size)
            text(0.6,0.4,['{ }$\alpha = ',num2str(alpha_vec(4)),'$'],'Units','normalized',...
                'Color',colors(4,:),'Interpreter','LaTeX','FontSize',font_size)
        end
        drawnow
    end

    if isequal(case_name,'2')
        figure(3); % reaction rate
        plot(r,k(r),linestyle{i},'Color',colors(i,:),'LineWidth',line_width)
        hold on
        if i == 1
            plot(r,kmin*ones(size(r)),'--','Color',[0.7,0.7,0.7],'LineWidth',line_width)
            plot(r,kmax*ones(size(r)),'--','Color',[0.7,0.7,0.7],'LineWidth',line_width)
        end
        ylim([kmin-0.05*(kmax-kmin),kmax+0.05*(kmax-kmin)])
        xlabel('$r$','Interpreter','LaTeX')
        set(gca,'Fontsize',font_size,'FontName','Times','Color',background_color,...
            'YTick',[kmin,kavg,kmax],'Layer','top','TickLabelInterpreter','LaTeX')
        labely = ylabel('$k(r)$','Interpreter','LaTeX');
        labely.Position(1) = -0.051*R; 
        labely.Position(2) = (kmin+kmax)/2;
        if i == 1
            text(-0.17,-0.17,'(b)','Units','normalized','Color','k','FontSize',font_size+12)
        end        
        set(gca,'YTickLabel',{'$k_{\rm min}$','$k_{\rm avg}$','$k_{\rm max}$'})
        set(gca,'XTickLabel',{'0','$R/2$','$R$'}')
        if i == 4
            text(0.6,0.7,['{ }$\alpha = ',num2str(alpha_vec(1)),'$'],'Units','normalized',...
                'Color',colors(1,:),'Interpreter','LaTeX','FontSize',font_size)
            text(0.6,0.6,['{ }$\alpha = ',num2str(alpha_vec(2)),'$'],'Units','normalized',...
                'Color',colors(2,:),'Interpreter','LaTeX','FontSize',font_size)
            text(0.6,0.5,['{ }$\alpha = ',num2str(alpha_vec(3)),'$'],'Units','normalized',...
                'Color',colors(3,:),'Interpreter','LaTeX','FontSize',font_size)
            text(0.6,0.4,['{ }$\alpha = ',num2str(alpha_vec(4)),'$'],'Units','normalized',...
                'Color',colors(4,:),'Interpreter','LaTeX','FontSize',font_size)
        end
        drawnow
    end

    if (i == 1 || i == 3) && cap_plots
        
        [X,Y,Z] = sphere(50);
        Xp = X(1:26,:);
        Yp = Y(1:26,:);
        Zp = Z(1:26,:);
        %cmap = cmap(end:-1:1,:);

        for j = 1:length(tsh)

            fh = figure;
            %set(gcf,'Renderer','Painters');
            for k = 2:length(rh)
                mapc = ceil(csh(k,j)*Nc);
                if rh(k) < 1
                    color_surf = cmap(mapc,:);
                else
                    color_surf = 0.7*ones(3,1);
                end
                surf(rh(k)*Xp,rh(k)*Yp,rh(k)*Zp,'EdgeColor','none','FaceColor',color_surf)%,'FaceAlpha',0.1)
                hold on
                surf(rh(k)*Xp,rh(k)*Zp,rh(k)*Yp,'EdgeColor','none','FaceColor',color_surf)
                %surf(k/N*Zp,k/N*Xp,k/N*Yp,'EdgeColor','none','FaceColor',cmap(k,:))
                axis equal
                view([120,21])
                xlim([-1,1])
                ylim([-1,1])
                zlim([-1,1])
                set(gca,'Xtick',[],'Ytick',[],'Ztick',[])
                text(2.9,1.5,['$\hat{t}=',num2str(tsh(j),'%g'),'$'],'FontSize',25,'FontName','Times','Interpreter','LaTeX')
                axis off
            end

            drawnow
            exportgraphics(fh,[fig_path,'conc_surf',case_name,num2str(i),num2str(j),'.pdf'],'Resolution',300)

        end

    end

end
fig2 = figure(2);
exportgraphics(fig2,[fig_path,'mass',case_name,'.pdf'])
fig3 = figure(3);
if isequal(case_name,'1')
    exportgraphics(fig3,[fig_path,'Dfunc.pdf'])
elseif isequal(case_name,'2')
    exportgraphics(fig3,[fig_path,'kfunc.pdf'])
end
