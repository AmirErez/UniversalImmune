clear all

dd2 = '../dat/';
fh = figure(1); clf
fh.Position = [360   278   560   600];
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
    
    subplot(3,2,r); hold on
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
            'fontsize',fs,'interpreter','latex'); 
        cb.Position = [0.42, 0.74, 0.03, 0.1];
    end
    
    if r == 2
        h_ = plot(10,10,'w',10,10,'w',10,10,'r-','linewidth',lw);
        leg = legend(h_,{'MEKi','SRCi','Hill'},'location','se',...
            'fontsize',fs2,'interpreter','latex');
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

delta_theta = 0.05;
[allSSEs, I1s, npoints, fitresult] = fitCv(delta_theta, D);

% Now estimate most likely
subplot(3,2,5);
ax = gca;


h = plot( fitresult, I1s, allSSEs);
% h = plot( fitresult, I1s', allSSEs, excludedPoints );
hold on

% fn = find((I1s>=0.2).*(I1s<=0.8));
% theory = fitresult.A + ((I1s(fn) - fitresult.x0).^2)/(2*fitresult.sigma^2);
% stddata = std(allSSEs(fn) - theory');

threshold = fitresult.A*(1+1/npoints);
xlim([0.25,0.75]);
plot([0.25,0.75], threshold*[1,1],'--k', 'DisplayName', 'Confidence');
dx = sqrt(fitresult.A*2*fitresult.sigma^2/npoints);
% halflen = round(length(allSSEs)/2);
% x0minus = interp1(allSSEs(1:halflen), I1s(1:halflen), threshold);
% x0plus = interp1(allSSEs(halflen:end), I1s(halflen:end), threshold);

xlabel('Scaling factor, $I_1$',...
    'fontsize',fs,'interpreter','latex')
ylabel('SSE',...
    'fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2)
text(0.26, 46700,'(e)','fontsize',fs3);
%t = text(0.37, 56000,{['$I_1=' num2str(fitresult.x0,2) '\pm' num2str(dx,2),'$']}, ...
%                      'fontsize',fs2, 'Interpreter', 'Latex');
box on
l = legend;
l.Interpreter = 'Latex';
l.Location = 'best';

% ----------- Now fits with different delta_theta
subplot(3,2,6);
delta_thetas = [0.03:0.01:0.1];
dIbars = zeros(size(delta_thetas));
Ibars = zeros(size(delta_thetas));
fitsigma = zeros(size(delta_thetas));
all_npoints = zeros(size(delta_thetas));
for ii=1:length(delta_thetas)
    disp(['delta_theta = ' num2str(delta_thetas(ii))]);
    [allSSEs, I1s, npoints, fitresult] = fitCv(delta_thetas(ii), D);
    Ibars(ii) = fitresult.x0;
    dIbars(ii) = sqrt(fitresult.A*2*fitresult.sigma^2/npoints);
    fitsigma(ii) = fitresult.sigma;
    all_npoints(ii) = npoints;
end
errorbar(delta_thetas, Ibars, dIbars);
xlabel('$\Delta \theta$', 'FontSize', fs2, 'Interpreter', 'latex');
ylabel('$I_1$', 'FontSize', fs2, 'Interpreter', 'latex');
ylim([0,0.8]);
t = text(0.005, 0.0952,'(f)','fontsize',fs3);

print(gcf,'-depsc',[fd 'fig4.eps'])

% Now plot supplement for npoints and fitsigma and vs. Delta theta
fh = figure();
yyaxis left
plot(delta_thetas, all_npoints, '.-');
ylabel('\# Experiment data points', 'FontSize', fs2, 'Interpreter', 'latex');
yyaxis right
plot(delta_thetas, 1./(fitsigma.^2), '.-');
ylabel('SSE curvature', 'FontSize', fs2, 'Interpreter', 'latex');
xlabel('$\Delta \theta$', 'FontSize', fs2, 'Interpreter', 'latex');
xlim([0.03,0.1]);
print(gcf,'-depsc',[fd 'fig4supp.eps'])