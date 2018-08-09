clear all

dd = '../dat/';
fh = figure(1); clf
fd = './';
fs = 16; fs2 = 12; fs3 = 20;
lw = 1.5; lw2 = 1;
gr = [0 .5 0];
pu = [.5 0 .5];

% free parameters
W = 25; % SG window (must be odd)
L = 100; % number of l bins

% free parameter used only for heat capacity
I1 = .1; % intensity of one molecule

% fixed parameter (minimum needed for derivatives)
J = 4; % SG polynomial degree

% cosmetic parameter (for plotting only)
NI = 500; % number of I bins

% Example experiment: MEKi (PD325901), intermediate dose

% loop over all experiments and drugs
expt = 'OT1Sig_20140919';
load([dd expt '.mat'])
v = 1;
i = 15;

% PDF of I
Ivals = expdata.RawData{v,i};
S = length(Ivals); % number of samples
I = linspace(min(Ivals),max(Ivals),NI)';
dI = I(2)-I(1);
ct = hist(Ivals,I)';
pI = ct/S/dI;

% PDF of n
nvals = Ivals/I1;
N = ceil(max(nvals));
n = (0:N)';
ct = hist(nvals,n)';
pn = ct/S;

% PDF of l
lvals = log(Ivals);
l = linspace(min(lvals),max(lvals),L)';
dl = l(2)-l(1);
ct = hist(lvals,l)';
q = ct/S/dl;

% SG filter and derivatives of Q
[B,G] = sgolay(J,W);
Y = (W+1)/2;
Q = -l+log(q);
for j = Y:L-Y
    q0(j-Y+1,1) = G(:,1)'*q(j-Y+1:j+Y-1);
    phi0(j-Y+1,1) = G(:,2)'*Q(j-Y+1:j+Y-1)/dl;
    phi1(j-Y+1,1) = 2*G(:,3)'*Q(j-Y+1:j+Y-1)/dl^2;
    phi2(j-Y+1,1) = 6*G(:,4)'*Q(j-Y+1:j+Y-1)/dl^3;
    phi3(j-Y+1,1) = 24*G(:,5)'*Q(j-Y+1:j+Y-1)/dl^4;
end
l0 = l(Y:L-Y);

% maxima (where phi=0); interpolate
j = find(diff(sign(phi0)) == -2);
lstar = (l0(j).*phi0(j+1)-l0(j+1).*phi0(j)) ...
    ./(phi0(j+1)-phi0(j));
nstar = exp(lstar)/I1;

% l_c: max of phi'-phi
phi10 = phi1 - phi0;
j = find(~isinf(phi10));
[ig,k] = max(phi10(j));
jc = j(k);
lc = l0(jc);

% theta, h
Ic = exp(lc)
theta = -2*phi1(jc)/(phi1(jc)-phi3(jc))
h = 2*phi0(jc)/(phi1(jc)-phi3(jc))


% ---
% A: q(l)
% ---

subplot(2,2,1); hold on
xl = [0 11];
yl = [0 .4];
hA1 = bar(l,q,'w','linewidth',lw2);
hA2 = plot(l0,q0,'r-','linewidth',lw);
xlim(xl)
ylim(yl)
xlabel('Log intensity, $\ell$','fontsize',fs,'interpreter','latex')
ylabel('Prob.\ density, $q(\ell)$',...
    'fontsize',fs,'interpreter','latex')
legend([hA1 hA2],{'Raw data','Filtered'},...
    'fontsize',fs2,'interpreter','latex',...
    'location','nw');
set(gca,'fontsize',fs2)
text(9,.34,'(a)','fontsize',fs3)
box on

% ---
% B: phi(l)
% ---

subplot(2,2,2);
dx = .2;
xl = [3 8];
yl = [-2.25 .75];
hB = plot(xl,[0 0],'k-',l0,phi0,'b-');
set(hB(1),'linewidth',lw2)
set(hB(2:end),'linewidth',lw)
xlim(xl)
ylim(yl)
xlabel('Log intensity, $\ell$','fontsize',fs,'interpreter','latex')
ylabel('$\phi(\ell)$','fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2)
text(3.2,-1.85,'(b)','fontsize',fs3)


% ---
% C: phi' - phi
% ---

subplot(2,2,3);
dx = .2;
xl = [3 8];
yl = [-2 4];
hC = plot(l0,phi1-phi0,'-',[lc lc],yl,'--');
set(hC,'linewidth',lw)
set(hC,'color',gr);
xlim(xl)
ylim(yl)
xlabel('Log intensity, $\ell$','fontsize',fs,'interpreter','latex')
ylabel('$\phi''-\phi$','fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2)
text(3.2,3.2,'(c)','fontsize',fs3)
text(lc+dx,.4,'$\ell_c$',...
    'fontsize',fs,'interpreter','latex','color',gr)


% ---
% D: p(I)
% ---

subplot(2,2,4)
dx = 170;
xl = [-100 4e3];
yl = [0 2.3e-3];
hD = plot(I,pI,'-',Ic*[1 1],yl,'--');
set(hD,'linewidth',lw)
set(hD(1),'color',pu)
set(hD(end),'color',gr)
xlim(xl)
ylim(yl)
xlabel('Intensity, $I$ (a.u.)',...
    'fontsize',fs,'interpreter','latex')
ylabel('Prob.\ density, $p(I)$',...
    'fontsize',fs,'interpreter','latex')
% legend(hD(end-1:end),...
%     {'$e^{\ell_*}$, from (b)','$e^{\ell_c}$, from (c)'},...
%     'fontsize',fs2,'interpreter','latex',...
%     'location','e');
set(gca,'fontsize',fs2)
text(Ic+dx,.95e-3,'$I_c$',...
    'fontsize',fs,'interpreter','latex','color',gr)
text(2e3,.95e-3,['$\theta = ' num2str(round(theta*100)/100) '$'],...
    'fontsize',fs,'interpreter','latex','color',pu)
text(2e3,.6e-3,['$h = ' num2str(round(h*100)/100) '$'],...
    'fontsize',fs,'interpreter','latex','color',pu)
text(3.3e3,2e-3,'(d)','fontsize',fs3)


print(gcf,'-depsc',[fd 'fig7.eps'])

