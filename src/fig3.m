clear all

dd2 = '../dat/';
fd = './';
fs = 16; fs2 = 12; fs3 = 20;
lw = 1.5; lw2 = 1;
ms = 8;

load([dd2 'fig3_data1.mat'])

% ---
% PDFs
% ---

fh1 = figure(1); clf
coA = [linspace(0,0,D(1))' linspace(0,1,D(1))' linspace(0,1,D(1))'];
coB = [linspace(0,1,D(2))' linspace(0,0,D(2))' linspace(0,1,D(2))'];

% MEKi
% plot PDFs of I
axA = subplot(2,2,1);
xl = [0 6e3];
yl = [0 1.6e-3];
hA = plot(Is{1}(:,end:-1:1),pIs{1}(:,end:-1:1),'-');
set(hA,'linewidth',lw)
xlim(xl)
ylim(yl)
for j = 1:D(1)
    set(hA(j),'color',coA(j,:));
end
xlabel('Intensity, $I$ (a.u.)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Prob.\ density, $p(I)$',...
    'fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2)
text(6e2,1.4e-3,'(a)','fontsize',fs3)

% plot PDFs of log I (inset)
a1 = get(fh1,'currentaxes');
a2 = axes('pos',[.25 .715 .2 .2]);
xl = [5e0 3e4];
yl = [0 .59];
hAi = semilogx(exp(ls{1}(:,end:-1:1)),qs{1}(:,end:-1:1),'-');
xlim(xl)
ylim(yl)
set(hAi,'linewidth',lw)
for j = 1:D(1)
    set(hAi(j),'color',coA(j,:));
end
xlabel('$I$','fontsize',fs,'interpreter','latex')
ylabel('$p$','fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2,'xtick',10.^(1:4))

% SRCi
% plot PDFs of I
axB = subplot(2,2,2);
xl = [0 6e3];
yl = [0 1.6e-3];
hB = plot(Is{2}(:,end:-1:1),pIs{2}(:,end:-1:1),'-');
set(hB,'linewidth',lw)
xlim(xl)
ylim(yl)
for j = 1:D(2)
    set(hB(j),'color',coB(j,:));
end
xlabel('Intensity, $I$ (a.u.)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Prob.\ density, $p(I)$',...
    'fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2)
text(6e2,1.4e-3,'(b)','fontsize',fs3)

% plot PDFs of log I (inset)
a1 = get(fh1,'currentaxes');
a2 = axes('pos',[.69 .715 .2 .2]);
xl = [5e0 3e4];
yl = [0 .59];
hBi = semilogx(exp(ls{2}(:,end:-1:1)),qs{2}(:,end:-1:1),'-');
xlim(xl)
ylim(yl)
set(hBi,'linewidth',lw)
for j = 1:D(2)
    set(hBi(j),'color',coB(j,:));
end
xlabel('$I$','fontsize',fs,'interpreter','latex')
ylabel('$p$','fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2,'xtick',10.^(1:4))


% ---
% "forcing" curves phi
% ---

% MEKi
axC = subplot(2,2,3);
xl = [2e1 3.5e3];
yl = [-3 1.5];
hC = semilogx(xl,[0 0],'k-',...
    exp(l0s{1}(:,end:-1:1)),phi0s{1}(:,end:-1:1),'-');
xlim(xl)
ylim(yl)
set(hC(1),'linewidth',lw2)
set(hC(2:end),'linewidth',lw)
for j = 1:D(1)
    set(hC(j+1),'color',coA(j,:));
end
xlabel('Intensity, $I$ (a.u.)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Forcing, $f(I) - I$',...
    'fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2,'xtick',10.^[2 3])
text(2.5e1,1,'(c)','fontsize',fs3)

% SRCi
axD = subplot(2,2,4);
xl = [2e1 4.2e3];
yl = [-4 2];
hD = semilogx(xl,[0 0],'k-',...
    exp(l0s{2}(:,end:-1:1)),phi0s{2}(:,end:-1:1),'-');
xlim(xl)
ylim(yl)
set(hD(1),'linewidth',lw2)
set(hD(2:end),'linewidth',lw)
for j = 1:D(2)
    set(hD(j+1),'color',coB(j,:));
end
xlabel('Intensity, $I$ (a.u.)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Forcing, $f(I) - I$',...
    'fontsize',fs,'interpreter','latex')
% set(get(cbB,'ylabel'),'string','[SRCi] (nM)',...
%     'fontsize',fs,'interpreter','latex')
%text(5.5e3,-5.1,{'[SRCi]','\,\,(nM)'},'fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2,'xtick',10.^[2 3])
text(2.7e1,1.3,'(d)','fontsize',fs3)

print(gcf,'-depsc',[fd 'fig3-1.eps'])


% ---
% theta, h (with error bars from filter windows), Ic (size)
% ---

load([dd2 'fig3_data2.mat'])
d = logspace(log10(1000),log10(1000/1.5^23),24);
Ic1 = 100;
Ic2 = 600;

figure(2); clf

% MEK colorbar
axA = subplot(2,2,1);
colormap(axA,coA)
cbC = colorbar;
log10maxdose = log10(1000);
log10mindose = log10(1000/1.5^23);
ticks = ((-1:3)-log10mindose)/(log10maxdose-log10mindose);
set(cbC,'location','north','xaxislocation','bottom',...
    'position',[.18 .75 .2 .04],...
    'ticks',ticks,'ticklabels',{'10^{-1}','10^0','10^1','10^2','10^3'},...
    'fontsize',fs2)
set(get(cbC,'title'),'string','[MEKi] (nM)',...
    'fontsize',fs,'interpreter','latex')

% SRC colorbar
axB = subplot(2,2,2);
colormap(axB,coB)
cbD = colorbar;
log10maxdose = log10(1000);
log10mindose = log10(1000/1.5^23);
ticks = ((-1:3)-log10mindose)/(log10maxdose-log10mindose);
set(cbD,'location','north','xaxislocation','bottom',...
    'position',[.63 .75 .2 .04],...
    'ticks',ticks,'ticklabels',{'10^{-1}','10^0','10^1','10^2','10^3'},...
    'fontsize',fs2)
set(get(cbD,'title'),'string','[SRCi] (nM)',...
    'fontsize',fs,'interpreter','latex')

% I_c, theta, h
subplot(2,8,9:10); hold on
xl = [6e-1 9e2];
yl = [0 650];
plot(xl,[0 0],'k-','linewidth',lw2)
for i = 1:D(1)
    plot(d(i)*[1 1],Icbars{1}(i)+dIcs{1}(i)*[-1 1],'k-',...
        'linewidth',lw2)
end
for i = 1:D(2)
    plot(d(i)*[1 1],Icbars{2}(i)+dIcs{2}(i)*[-1 1],'k-',...
        'linewidth',lw2)
end
for i = 1:D(1)
    plot(d(D(1)-i+1),Icbars{1}(D(1)-i+1),'o',...
        'color',coA(i,:),'markersize',ms,'markerfacecolor',coA(i,:));
end
for i = 1:D(2)
    plot(d(D(2)-i+1),Icbars{2}(D(2)-i+1),'s',...
        'color',coB(i,:),'markersize',ms,'markerfacecolor',coB(i,:));
end
xlim(xl)
ylim(yl)
xlabel('Dose (nM)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Intensity, $I_c$',...
    'fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2,'xtick',10.^(-1:3),'xscale','log')
text(1.2e0,75,'(e)','fontsize',fs3)
box on

subplot(2,8,12:13); hold on
yl = [-1 .5];
%plot(xl,[0 0],'k-','linewidth',lw2)
for i = 1:D(1)
    plot(d(i)*[1 1],thetabars{1}(i)+dthetas{1}(i)*[-1 1],'k-',...
        'linewidth',lw2)
end
for i = 1:D(2)
    plot(d(i)*[1 1],thetabars{2}(i)+dthetas{2}(i)*[-1 1],'k-',...
        'linewidth',lw2)
end
for i = 1:D(1)
    plot(d(D(1)-i+1),thetabars{1}(D(1)-i+1),'o',...
        'color',coA(i,:),'markersize',ms,'markerfacecolor',coA(i,:));
end
for i = 1:D(2)
    plot(d(D(2)-i+1),thetabars{2}(D(2)-i+1),'s',...
        'color',coB(i,:),'markersize',ms,'markerfacecolor',coB(i,:));
end
xlim(xl)
ylim(yl)
xlabel('Dose (nM)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Reduced temp., $\theta$',...
    'fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2,'xtick',10.^(-1:3),'xscale','log')
text(1.2e0,.3,'(f)','fontsize',fs3)
box on

subplot(2,8,15:16); hold on
yl = [-1.3 -.1];
%plot(xl,[0 0],'k-','linewidth',lw2)
for i = 1:D(1)
    plot(d(i)*[1 1],hbars{1}(i)+dhs{1}(i)*[-1 1],'k-',...
        'linewidth',lw2)
end
for i = 1:D(2)
    plot(d(i)*[1 1],hbars{2}(i)+dhs{2}(i)*[-1 1],'k-',...
        'linewidth',lw2)
end
for i = 1:D(1)
    plot(d(D(1)-i+1),hbars{1}(D(1)-i+1),'o',...
        'color',coA(i,:),'markersize',ms,'markerfacecolor',coA(i,:));
end
for i = 1:D(2)
    plot(d(D(2)-i+1),hbars{2}(D(2)-i+1),'s',...
        'color',coB(i,:),'markersize',ms,'markerfacecolor',coB(i,:));
end
xlim(xl)
ylim(yl)
xlabel('Dose (nM)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Field, $h$',...
    'fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2,'xtick',10.^(-1:3),'xscale','log')
text(1.2e0,-1.15,'(g)','fontsize',fs3)
box on

print(gcf,'-depsc',[fd 'fig3-2.eps'])


