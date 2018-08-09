clear all

dd2 = '../dat/';
fh = figure(1); clf
fd = './';
fs = 16; fs2 = 12; fs3 = 20;
lw = 1.5;
ms1 = 12; ms2 = 5;
br = [152 102 51]/255;
ch = 'abcd';

% MEK: ^ PD325901, v AZD6244, > Trametinib
% SRC: o Dasatinib, s Bosutinib, d PP2
sh = '^o^ov>^sd';

% free parameter used only for heat capacity
I1s = [10 1 .1 .01]; % intensity of one molecule

Z = 1e3;
co = colormap(parula(Z));
hmin = -.6;
hmax = 0;
hmap = linspace(hmin,hmax,Z);

for r = 1:length(I1s)
    
    load([dd2 'fig4_data1_I' num2str(r) '.mat'])
    
    subplot(2,2,r); hold on
    for w = 1:length(D)
        ms = linspace(ms1,ms2,D(w));
        for i = 1:D(w)
            [ig,j] = min(abs(hs{w}(i)-hmap));
            plot(thetas{w}(i),Cs{w}(i),sh(w),'markersize',ms(i),...
                'markerfacecolor',co(j,:),'markeredgecolor','k');
        end
    end
    [ig,j] = min(abs(h-hmap));
    hD = plot(tH,CH,'r-','linewidth',lw);
    
    if r == 1
        cb = colorbar;
        nt = 4;
        set(cb,'ticks',linspace(0,1,nt),...
            'ticklabels',linspace(hmin,hmax,nt),...
            'yaxislocation','left','fontsize',fs2,...
            'position',[.42 .62 .03 .18])
        set(get(cb,'ylabel'),'string','Field, $h$',...
            'fontsize',fs,'interpreter','latex')
    end
    
    if r == 2
        h_ = plot(10,10,'w',10,10,'w',10,10,'r-','linewidth',lw);
        legend(h_,{'MEKi','SRCi','Hill'},'location','se',...
            'fontsize',fs2,'interpreter','latex')
    end
    
    xlim(xl)
    ylim(yl)
    xlabel('Reduced temperature, $\theta$',...
        'fontsize',fs,'interpreter','latex')
    ylabel('Heat capacity, $C/k_{\rm B}$',...
        'fontsize',fs,'interpreter','latex')
    title(['$I_1 = ' num2str(I1) '$'],...
        'fontsize',fs,'interpreter','latex')
    set(gca,'fontsize',fs2)
    text(-.9,ys(r),['(' ch(r) ')'],'fontsize',fs3)
    box on
end

print(gcf,'-depsc',[fd 'fig4.eps'])

