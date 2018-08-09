clear all

fh = figure(1); clf
fd = './';
fs = 16; fs2 = 12; fs3 = 14; fs4 = 20;
lw = 1.5; lw2 = 1;
ms = 6; mhs = 3;
gr = [0 .5 0];
pu = [.5 0 .5];
or = [1 .5 0];
ye = [0.9290 0.6940 0.1250];
br = [153 102 51]/256;
lb = [0.3010 0.7450 0.9330];
ma = [0.6350 0.0780 0.1840];


% ---
% A: Cartoon
% ---

% ---
% B: forcing function (Schlogl, Hill)
% ---

nc = 100;
N = 6*nc;
n = (0:N)';

% Hill (a varies with s to keep h = 0)
Z = 100;
H = 3
K = nc*((H+1)/(H-1))^(1/H)
sc = 4*H*nc/(H^2-1);
ss = linspace(.3*sc,1.7*sc,Z);
ac = nc - (H-1)*sc/2/H;
as = nc - (H-1)*ss/2/H;
s = ss(Z) % farthest into bistable regime
a = as(Z)
fH = a+s*n.^H./(n.^H+K^H);
j = find(diff(sign(fH-n)) == -2);
nsH = n(j);
j1 = find(diff(sign(fH-n)) == 2);
j2 = max(find(diff(sign(fH-n)) == 1));
nuH = n([j1 j2]);

% Schlogl (same parameters as Hill)
fS = a*K^2./((n-1).*(n-2)+K^2)+s*(n-1).*(n-2)./((n-1).*(n-2)+K^2);
j = find(diff(sign(fS-n)) == -2);
nsS = n(j);
j1 = find(diff(sign(fS-n)) == 2);
j2 = max(find(diff(sign(fS-n)) == 1));
nuS = n([j1 j2]);

% plot feedback functions (f_n - n)
subplot(2,2,2)
xl = [0 250];
x1 = 100;
dx = 12;
y1 = -16;
y2 = -22;
hB = plot(xl,[0 0],'k-',n,fH-n,'r-',n,fS-n,'b-',...
    nsH,[0 0],'ro',nsS,[0 0],'bo',nuH,0,'ro',nuS,0,'bo',...
    x1,y1,'ko',x1,y2,'ko');
set(hB,'linewidth',lw2)
set(hB(2:3),'linewidth',lw)
set(hB(4:end),'markersize',ms)
set(hB(4),'markerfacecolor','r')
set(hB(5),'markerfacecolor','b')
set(hB([6 7 9]),'markerfacecolor','w')
set(hB(8),'markerfacecolor','k')
xlim(xl)
ylim([-30 40])
xlabel('Number of molecules, $n$',...
    'fontsize',fs,'interpreter','latex')
ylabel('Forcing, $f_n-n$','fontsize',fs,'interpreter','latex')
legend(hB([3 2]),{'Schl\"ogl','Hill'},...
    'location','nw','fontsize',fs2,'interpreter','latex')
text(x1+dx,y1,'stable','fontsize',fs2,'interpreter','latex')
text(x1+dx,y2,'unstable','fontsize',fs2,'interpreter','latex')
set(gca,'fontsize',fs2)
text(210,32,'(b)','fontsize',fs4)


% ---
% D: distributions for varying h (and theta legend)
% ---

theta = 0;
% h < 0
h = -.03;
s = nc*16*H/((H^2-1)*((H^2-1)*theta+4));
a = nc*(H-1)*((H+1)^2*(theta+h)+4)/((H+1)*((H^2-1)*theta+4));
f = a+s*n.^H./(n.^H+K^H);
c = [1; cumprod(f(2:end)./n(2:end))];
p1 = c/sum(c);
lstrD{1} = ['$h = ' num2str(h) '$'];

% h = 0
h = 0;
s = nc*16*H/((H^2-1)*((H^2-1)*theta+4));
a = nc*(H-1)*((H+1)^2*(theta+h)+4)/((H+1)*((H^2-1)*theta+4));
f = a+s*n.^H./(n.^H+K^H);
c = [1; cumprod(f(2:end)./n(2:end))];
p2 = c/sum(c);
lstrD{2} = ['$h = ' num2str(h) '$'];

% h > 0
h = .03;
s = nc*16*H/((H^2-1)*((H^2-1)*theta+4));
a = nc*(H-1)*((H+1)^2*(theta+h)+4)/((H+1)*((H^2-1)*theta+4));
f = a+s*n.^H./(n.^H+K^H);
c = [1; cumprod(f(2:end)./n(2:end))];
p3 = c/sum(c);
lstrD{3} = ['$h = ' num2str(h) '$'];

% plot distributions
subplot(2,2,4)
hD = plot([nc nc],[0 .04],'k--',n,p1,'-',n,p2,'-',n,p3,'-');
set(hD,'linewidth',lw)
set(hD(2),'color',lb)
set(hD(3),'color',gr)
set(hD(4),'color',br)
xlim([0 260])
ylim([0 .04])
xlabel('Number of molecules, $n$',...
    'fontsize',fs,'interpreter','latex')
ylabel('Probability, $p_n$',...
    'fontsize',fs,'interpreter','latex')
text(107,.019,'$n_c$','fontsize',fs,'interpreter','latex')
legend(hD(2:4),lstrD,...
    'location',[.81 .375 eps eps],...
    'fontsize',fs2,'interpreter','latex')
set(gca,'fontsize',fs2,'xtick',0:50:250)
text(220,.005,'(d)','fontsize',fs4)


% ---
% C: distributions and pitchfork (Hill, a varies with s to keep h = 0)
% ---

% unimodal
s1 = 100;
s = s1;
a = nc - (H-1)*s/2/H;
f = a+s*n.^H./(n.^H+K^H);
c = [1; cumprod(f(2:end)./n(2:end))];
p1 = c/sum(c);
f1 = (H^2-1)*s/4/H/nc;
f3 = -(H^2-1)^2*s/8/H/nc^3;
theta = -2*(1-f1)/f3/nc^2;
lstrC{1} = ['$\theta = ' num2str(theta) '$'];

% critical
s = sc;
a = ac;
f = a+s*n.^H./(n.^H+K^H);
c = [1; cumprod(f(2:end)./n(2:end))];
p2 = c/sum(c);
f1 = (H^2-1)*s/4/H/nc;
f3 = -(H^2-1)^2*s/8/H/nc^3;
theta = -2*(1-f1)/f3/nc^2;
lstrC{2} = ['$\theta = ' num2str(theta) '$'];

% bimodal
s2 = 200;
s = s2;
a = nc - (H-1)*s/2/H;
f = a+s*n.^H./(n.^H+K^H);
c = [1; cumprod(f(2:end)./n(2:end))];
p3 = c/sum(c);
f1 = (H^2-1)*s/4/H/nc;
f3 = -(H^2-1)^2*s/8/H/nc^3;
theta = -2*(1-f1)/f3/nc^2;
lstrC{3} = ['$\theta = ' num2str(theta) '$'];

% plot distributions
subplot(2,2,3)
hC = plot(n,p1,'-',n,p2,'-',n,p3,'-');
set(hC,'linewidth',lw)
set(hC(1),'color',or)
set(hC(2),'color',gr)
set(hC(3),'color',pu)
xlim([0 260])
ylim([0 .04])
xlabel('Number of molecules, $n$',...
    'fontsize',fs,'interpreter','latex')
ylabel('Probability, $p_n\qquad\qquad$',...
    'fontsize',fs,'interpreter','latex')
legend(hC(end:-1:1),lstrC(end:-1:1),...
    'location',[.14 .41 eps eps],...
    'fontsize',fs2,'interpreter','latex')
set(gca,'fontsize',fs2,'xtick',0:50:250)
text(220,.005,'(c)','fontsize',fs4)

% Hill pitchfork
n_ = (0:.1:N/2)';
n1 = []; n2 = []; n3 = [];
for i = 1:Z
    a = as(i);
    s = ss(i);
    f = a+s*n_.^H./(n_.^H+K^H);
    j = find(diff(sign(f-n_)) < 0); % sign + -> - (stable)
    k = find(diff(sign(f-n_)) > 0); % sign - -> + (unstable)
    n1(i) = min(n_(j));
    n2(i) = max(n_(j));
    if isempty(k)
        n3(i) = n1(i);
    else
        n3(i) = min(n_(k));
    end
end

% plot pitchfork (inset)
a1 = get(fh,'currentaxes');
a2 = axes('pos',[.33 .265 .16*[.8 1.1]]);
hold on
hBi = plot(ss,n3,'r--',ss,n1,'r-',ss,n2,'r-');
set(hBi,'linewidth',lw)
p1 = [s1 0]; p2 = [s1 -40]; dp = p2-p1;
quiver(p1(1),p1(2),dp(1),dp(2),0,'maxheadsize',mhs,...
    'color',or,'linewidth',lw)
p1 = [sc 0]; p2 = [sc -40]; dp = p2-p1;
quiver(p1(1),p1(2),dp(1),dp(2),0,'maxheadsize',mhs,...
    'color',gr,'linewidth',lw)
p1 = [s2 0]; p2 = [s2 -40]; dp = p2-p1;
quiver(p1(1),p1(2),dp(1),dp(2),0,'maxheadsize',mhs,...
    'color',pu,'linewidth',lw)
xlim([75 225])
ylim([-50 225])
xlabel('$s$','fontsize',fs,'interpreter','latex')
ylabel('$n_*$','fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2)
box on



print(gcf,'-depsc',[fd 'fig1.eps'])

