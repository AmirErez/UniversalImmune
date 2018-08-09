clear all

fh = figure(1); clf
fd = './';
fs = 16; fs2 = 12; fs3 = 14; fs4 = 20;
lw = 1.5; lw2 = 1;
ms = 6;
dx = .44; dy = .475;

nc = 500;
N = 2.5*nc;
n = (0:N)';
Z1 = 4*10;
Z2 = 100;
Z3 = 20;
Z4 = 2000;

% ---
% A: beta (m vs theta at h = 0)
% ---

h = 0;
theta1 = logspace(-4,-1,Z1);
theta2 = linspace(-.2,.2,Z2);

% main panel: log in theta < 0
for i = 1:Z1
    theta = -theta1(i);
    
    % Schlogl
    x = 2*nc-3;
    K = sqrt(3*x^2+1)/2;
    s = (3*nc^3*(theta+h)+nc*x^2+x^3)/(3*nc^2*theta+x^2);
    a = ((3*x^2+1)*(3*nc^3*(theta+h)+nc*x^2+x^3)-4*x^5)...
        /(3*x^2+1)/(3*nc^2*theta+x^2);
    f = a*K^2./((n-1).*(n-2)+K^2)+s*(n-1).*(n-2)./((n-1).*(n-2)+K^2);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1S(i) = min(n(j));
    n2S(i) = max(n(j));
    
    % Hill
    H = 3;
    K = nc*((H+1)/(H-1))^(1/H);
    s = 16*H*nc/(H^2-1)/((H^2-1)*theta+4);
    a = nc*(H-1)*((H+1)^2*(theta+h)+4)/(H+1)/((H^2-1)*theta+4);
    f = a+s*n.^H./(n.^H+K^H);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1H(i) = min(n(j));
    n2H(i) = max(n(j));
end

% order parameter
m1S = (n1S-nc)/nc;
m2S = (n2S-nc)/nc;
m1H = (n1H-nc)/nc;
m2H = (n2H-nc)/nc;

% small m prediction
beta = 1/2;
m_ = (3*theta1).^beta;

% plot
j1 = 1:4:Z1;
j2 = 2:4:Z1;
j3 = 3:4:Z1;
j4 = 4:4:Z1;

subplot(2,2,1)
hA = loglog(theta1,m_,'k-',...
    theta1(j1),-m1S(j1),'bo',...
    theta1(j2),m2S(j2),'bs',...
    theta1(j3),-m1H(j3),'r^',...
    theta1(j4),m2H(j4),'rv');
set(hA,'linewidth',lw,'markersize',ms)
set(hA(3),'markerfacecolor','b')
set(hA(5),'markerfacecolor','r')
xlim([min(theta1) max(theta1)])
ylim([1e-2 7.3e-1])
xlabel('Neg.\ reduced temp., $-\theta$',...
    'fontsize',fs,'interpreter','latex')
ylabel('Abs.\ order param., $|m|$',...
    'fontsize',fs,'interpreter','latex')
legend(hA([2:end 1]),{'Schl\"ogl, $m<0$','Schl\"ogl, $m>0$',...
    'Hill, $m<0$','Hill, $m>0$','$\beta = 1/2$'},...
    'fontsize',fs2,'interpreter','latex',...
    'location',[.19 .89 eps eps])
set(gca,'fontsize',fs2)
text(4.5e-3,4.5e-1,'(a)','fontsize',fs4)


% inset: linear in theta
n1S = []; n2S = []; n1H = []; n2H = [];
for i = 1:Z2
    theta = theta2(i);
    
    % Schlogl
    x = 2*nc-3;
    K = sqrt(3*x^2+1)/2;
    s = (3*nc^3*(theta+h)+nc*x^2+x^3)/(3*nc^2*theta+x^2);
    a = ((3*x^2+1)*(3*nc^3*(theta+h)+nc*x^2+x^3)-4*x^5)...
        /(3*x^2+1)/(3*nc^2*theta+x^2);
    f = a*K^2./((n-1).*(n-2)+K^2)+s*(n-1).*(n-2)./((n-1).*(n-2)+K^2);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1S(i) = min(n(j));
    n2S(i) = max(n(j));
    
    % Hill
    H = 3;
    K = nc*((H+1)/(H-1))^(1/H);
    s = 16*H*nc/(H^2-1)/((H^2-1)*theta+4);
    a = nc*(H-1)*((H+1)^2*(theta+h)+4)/(H+1)/((H^2-1)*theta+4);
    f = a+s*n.^H./(n.^H+K^H);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1H(i) = min(n(j));
    n2H(i) = max(n(j));
end

% order parameter
m1S = (n1S-nc)/nc;
m2S = (n2S-nc)/nc;
m1H = (n1H-nc)/nc;
m2H = (n2H-nc)/nc;

% only plot one curve at m = 0 for t > 0
j = find(theta2 > 0);
m1S(j) = 0;
m1H(j) = 0;
j = find(theta2 <= 0);
theta2_ = theta2(j);
m2S_ = m2S(j);
m2H_ = m2H(j);

% plot inset
a1 = get(fh,'currentaxes');
a2 = axes('pos',[.34 .61 .13*[.8 1]]);
hAi = plot(theta2,m1S,'b-',theta2_,m2S_,'b-',...
    theta2,m1H,'r--',theta2_,m2H_,'r--');
set(hAi,'linewidth',lw)
xlim([min(theta2) max(theta2)])
ylim([-1 1.5])
xlabel('$\theta$','fontsize',fs,'interpreter','latex')
ylabel('$m$','fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2,'xaxislocation','top',...
    'xtick',[-.2 0 .2])


% ---
% B: gamma (chi vs theta at h = 0)
% ---

theta1 = logspace(-4,-1,Z1);
theta2 = linspace(-.2,.2,Z2);

% use adaptive range for slope estimation
dh1 = logspace(-6,-2,Z1);
dh2 = [linspace(1e-2,1e-4,Z2/2) linspace(1e-4,1e-2,Z2/2)];

% main panel: log in theta
% theta < 0
for i = 1:Z1
    theta = -theta1(i);
    
    % Schlogl
    x = 2*nc-3;
    K = sqrt(3*x^2+1)/2;
    
    % h = 0-
    h = -dh1(i)/2;
    s = (3*nc^3*(theta+h)+nc*x^2+x^3)/(3*nc^2*theta+x^2);
    a = ((3*x^2+1)*(3*nc^3*(theta+h)+nc*x^2+x^3)-4*x^5)...
        /(3*x^2+1)/(3*nc^2*theta+x^2);
    f = a*K^2./((n-1).*(n-2)+K^2)+s*(n-1).*(n-2)./((n-1).*(n-2)+K^2);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1nSa = min(n(j));
    n2nSa = max(n(j));
    
    % h = 0+
    h = dh1(i)/2;
    s = (3*nc^3*(theta+h)+nc*x^2+x^3)/(3*nc^2*theta+x^2);
    a = ((3*x^2+1)*(3*nc^3*(theta+h)+nc*x^2+x^3)-4*x^5)...
        /(3*x^2+1)/(3*nc^2*theta+x^2);
    f = a*K^2./((n-1).*(n-2)+K^2)+s*(n-1).*(n-2)./((n-1).*(n-2)+K^2);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1nSb = min(n(j));
    n2nSb = max(n(j));
    dn1nS(i) = n1nSb-n1nSa;
    dn2nS(i) = n2nSb-n2nSa;

    % Hill
    H = 3;
    K = nc*((H+1)/(H-1))^(1/H);
    s = 16*H*nc/(H^2-1)/((H^2-1)*theta+4);

    % h = 0-
    h = -dh1(i)/2;
    a = nc*(H-1)*((H+1)^2*(theta+h)+4)/(H+1)/((H^2-1)*theta+4);
    f = a+s*n.^H./(n.^H+K^H);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1nHa = min(n(j));
    n2nHa = max(n(j));

    % h = 0+
    h = dh1(i)/2;
    a = nc*(H-1)*((H+1)^2*(theta+h)+4)/(H+1)/((H^2-1)*theta+4);
    f = a+s*n.^H./(n.^H+K^H);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1nHb = min(n(j));
    n2nHb = max(n(j));
    dn1nH(i) = n1nHb-n1nHa;
    dn2nH(i) = n2nHb-n2nHa;
end

% theta > 0
for i = 1:Z1
    theta = theta1(i);
    
    % Schlogl
    x = 2*nc-3;
    K = sqrt(3*x^2+1)/2;
    
    % h = 0-
    h = -dh1(i)/2;
    s = (3*nc^3*(theta+h)+nc*x^2+x^3)/(3*nc^2*theta+x^2);
    a = ((3*x^2+1)*(3*nc^3*(theta+h)+nc*x^2+x^3)-4*x^5)...
        /(3*x^2+1)/(3*nc^2*theta+x^2);
    f = a*K^2./((n-1).*(n-2)+K^2)+s*(n-1).*(n-2)./((n-1).*(n-2)+K^2);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1pSa = min(n(j));
    n2pSa = max(n(j));
    
    % h = 0+
    h = dh1(i)/2;
    s = (3*nc^3*(theta+h)+nc*x^2+x^3)/(3*nc^2*theta+x^2);
    a = ((3*x^2+1)*(3*nc^3*(theta+h)+nc*x^2+x^3)-4*x^5)...
        /(3*x^2+1)/(3*nc^2*theta+x^2);
    f = a*K^2./((n-1).*(n-2)+K^2)+s*(n-1).*(n-2)./((n-1).*(n-2)+K^2);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1pSb = min(n(j));
    n2pSb = max(n(j));
    dn1pS(i) = n1pSb-n1pSa;
    dn2pS(i) = n2pSb-n2pSa;

    % Hill
    H = 3;
    K = nc*((H+1)/(H-1))^(1/H);
    s = 16*H*nc/(H^2-1)/((H^2-1)*theta+4);

    % h = 0-
    h = -dh1(i)/2;
    a = nc*(H-1)*((H+1)^2*(theta+h)+4)/(H+1)/((H^2-1)*theta+4);
    f = a+s*n.^H./(n.^H+K^H);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1pHa = min(n(j));
    n2pHa = max(n(j));

    % h = 0+
    h = dh1(i)/2;
    a = nc*(H-1)*((H+1)^2*(theta+h)+4)/(H+1)/((H^2-1)*theta+4);
    f = a+s*n.^H./(n.^H+K^H);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1pHb = min(n(j));
    n2pHb = max(n(j));
    dn1pH(i) = n1pHb-n1pHa;
    dn2pH(i) = n2pHb-n2pHa;
end

% susceptibility (derivative of order parameter)
chi1nS = dn1nS./dh1/nc;
chi2nS = dn2nS./dh1/nc;
chi1nH = dn1nH./dh1/nc;
chi2nH = dn2nH./dh1/nc;
chi1pS = dn1pS./dh1/nc;
chi2pS = dn2pS./dh1/nc;
chi1pH = dn1pH./dh1/nc;
chi2pH = dn2pH./dh1/nc;

% small m prediction
gamma_ = 1;
chi_n = (2*theta1).^(-gamma_);
chi_p = theta1.^(-gamma_);

% plot
j1 = 1:4:Z1;
j2 = 2:4:Z1;
j3 = 3:4:Z1;
j4 = 4:4:Z1;

subplot(2,2,2)
hB = loglog(theta1,chi_n,'k-',...
    theta1,chi_p,'k-',...
    theta1(j1),chi1nS(j1),'bo',...
    theta1(j2),chi1pS(j2),'bs',...
    theta1(j3),chi1nH(j3),'r^',...
    theta1(j4),chi1pH(j4),'rv');
set(hB,'linewidth',lw,'markersize',ms)
set(hB(4),'markerfacecolor','b')
set(hB(6),'markerfacecolor','r')
xlim([min(theta1) max(theta1)])
ylim([2e0 1e4])
xlabel('Abs. reduced temp., $|\theta|$',...
    'fontsize',fs,'interpreter','latex')
ylabel('Susceptibility, $\chi$',...
    'fontsize',fs,'interpreter','latex')
legend(hB([3:end 1]),{'Schl\"ogl, $\theta<0$',...
    'Schl\"ogl, $\theta>0$',...
    'Hill, $\theta<0$',...
    'Hill, $\theta>0$',...
    '$\gamma = 1$'},...
    'fontsize',fs2,'interpreter','latex',...
    'location',[.825 .41+dy eps eps])
set(gca,'fontsize',fs2)
text(4e-3,1.2e1,'$\chi$','rotation',90,...
    'fontsize',fs,'interpreter','latex')
text(7e-4,4e3,'(b)','fontsize',fs4)

% inset: linear in theta
for i = 1:Z2
    theta = theta2(i);
    
    % Schlogl
    x = 2*nc-3;
    K = sqrt(3*x^2+1)/2;
    
    % h = 0-
    h = -dh2(i)/2;
    s = (3*nc^3*(theta+h)+nc*x^2+x^3)/(3*nc^2*theta+x^2);
    a = ((3*x^2+1)*(3*nc^3*(theta+h)+nc*x^2+x^3)-4*x^5)...
        /(3*x^2+1)/(3*nc^2*theta+x^2);
    f = a*K^2./((n-1).*(n-2)+K^2)+s*(n-1).*(n-2)./((n-1).*(n-2)+K^2);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1Sa = min(n(j));
    n2Sa = max(n(j));
    
    % h = 0+
    h = dh2(i)/2;
    s = (3*nc^3*(theta+h)+nc*x^2+x^3)/(3*nc^2*theta+x^2);
    a = ((3*x^2+1)*(3*nc^3*(theta+h)+nc*x^2+x^3)-4*x^5)...
        /(3*x^2+1)/(3*nc^2*theta+x^2);
    f = a*K^2./((n-1).*(n-2)+K^2)+s*(n-1).*(n-2)./((n-1).*(n-2)+K^2);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1Sb = min(n(j));
    n2Sb = max(n(j));
    dn1S(i) = n1Sb-n1Sa;
    dn2S(i) = n2Sb-n2Sa;

    % Hill
    H = 3;
    K = nc*((H+1)/(H-1))^(1/H);
    s = 16*H*nc/(H^2-1)/((H^2-1)*theta+4);

    % h = 0-
    h = -dh2(i)/2;
    a = nc*(H-1)*((H+1)^2*(theta+h)+4)/(H+1)/((H^2-1)*theta+4);
    f = a+s*n.^H./(n.^H+K^H);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1Ha = min(n(j));
    n2Ha = max(n(j));

    % h = 0+
    h = dh2(i)/2;
    a = nc*(H-1)*((H+1)^2*(theta+h)+4)/(H+1)/((H^2-1)*theta+4);
    f = a+s*n.^H./(n.^H+K^H);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1Hb = min(n(j));
    n2Hb = max(n(j));
    dn1H(i) = n1Hb-n1Ha;
    dn2H(i) = n2Hb-n2Ha;

end

% susceptibility (derivative of order parameter)
chi1S = dn1S./dh2/nc;
chi2S = dn2S./dh2/nc;
chi1H = dn1H./dh2/nc;
chi2H = dn2H./dh2/nc;

% plot inset
a1 = get(fh,'currentaxes');
a2 = axes('pos',[.5925 .605 .13*[.8 1]]);
hBi = plot(theta2,chi1S,'b-',theta2,chi1H,'r--');
set(hBi,'linewidth',lw)
xlim([min(theta2) max(theta2)])
ylim([0 100])
xlabel('$\theta$','fontsize',fs,'interpreter','latex')
%ylabel('$\chi$','fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2,'xaxislocation','top',...
    'yaxislocation','right','xtick',[-.2 0 .2])


% ---
% C: delta (m vs h at theta = 0)
% ---

theta = 0;
h1 = logspace(-4,-1,Z1);
h2 = linspace(-.2,.2,Z2);

% main panel: log in h
% h < 0
n1S = []; n2S = []; n1H = []; n2H = [];
for i = 1:Z1
    h = -h1(i);
    
    % Schlogl
    x = 2*nc-3;
    K = sqrt(3*x^2+1)/2;
    s = (3*nc^3*(theta+h)+nc*x^2+x^3)/(3*nc^2*theta+x^2);
    a = ((3*x^2+1)*(3*nc^3*(theta+h)+nc*x^2+x^3)-4*x^5)...
        /(3*x^2+1)/(3*nc^2*theta+x^2);
    f = a*K^2./((n-1).*(n-2)+K^2)+s*(n-1).*(n-2)./((n-1).*(n-2)+K^2);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1S(i) = mean(n(j));
    
    % Hill
    H = 3;
    K = nc*((H+1)/(H-1))^(1/H);
    s = 16*H*nc/(H^2-1)/((H^2-1)*theta+4);
    a = nc*(H-1)*((H+1)^2*(theta+h)+4)/(H+1)/((H^2-1)*theta+4);
    f = a+s*n.^H./(n.^H+K^H);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n1H(i) = mean(n(j));
end

% h > 0
for i = 1:Z1
    h = h1(i);
    
    % Schlogl
    x = 2*nc-3;
    K = sqrt(3*x^2+1)/2;
    s = (3*nc^3*(theta+h)+nc*x^2+x^3)/(3*nc^2*theta+x^2);
    a = ((3*x^2+1)*(3*nc^3*(theta+h)+nc*x^2+x^3)-4*x^5)...
        /(3*x^2+1)/(3*nc^2*theta+x^2);
    f = a*K^2./((n-1).*(n-2)+K^2)+s*(n-1).*(n-2)./((n-1).*(n-2)+K^2);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n2S(i) = mean(n(j));
    
    % Hill
    H = 3;
    K = nc*((H+1)/(H-1))^(1/H);
    s = 16*H*nc/(H^2-1)/((H^2-1)*theta+4);
    a = nc*(H-1)*((H+1)^2*(theta+h)+4)/(H+1)/((H^2-1)*theta+4);
    f = a+s*n.^H./(n.^H+K^H);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    n2H(i) = mean(n(j));
end

% order parameter
m1S = (n1S-nc)/nc;
m2S = (n2S-nc)/nc;
m1H = (n1H-nc)/nc;
m2H = (n2H-nc)/nc;

% small m prediction
delta = 3;
m_ = (3*h1).^(1/delta);

% plot
j1 = 1:4:Z1;
j2 = 2:4:Z1;
j3 = 3:4:Z1;
j4 = 4:4:Z1;

subplot(2,4,6:7)
hC = loglog(h1,m_,'k-',...
    h1(j1),-m1S(j1),'bo',...
    h1(j2),m2S(j2),'bs',...
    h1(j3),-m1H(j3),'r^',...
    h1(j4),m2H(j4),'rv');
set(hC,'linewidth',lw,'markersize',ms)
set(hC(3),'markerfacecolor','b')
set(hC(5),'markerfacecolor','r')
xlim([min(h1) max(h1)])
ylim([3e-2 1])
xlabel('Absolute field, $|h|$',...
    'fontsize',fs,'interpreter','latex')
ylabel('Abs.\ order param., $|m|$',...
    'fontsize',fs,'interpreter','latex')
legend(hC([2:end 1]),{'Schl\"ogl, $h<0$','Schl\"ogl, $h>0$',...
    'Hill, $h<0$','Hill, $h>0$','$\delta = 3$'},...
    'fontsize',fs2,'interpreter','latex',...
    'location',[.39 .4 eps eps])
set(gca,'fontsize',fs2)
text(4e-3,6.5e-1,'(c)','fontsize',fs4)

% inset: linear in h
for i = 1:Z2
    h = h2(i);
    
    % Schlogl
    x = 2*nc-3;
    K = sqrt(3*x^2+1)/2;
    s = (3*nc^3*(theta+h)+nc*x^2+x^3)/(3*nc^2*theta+x^2);
    a = ((3*x^2+1)*(3*nc^3*(theta+h)+nc*x^2+x^3)-4*x^5)...
        /(3*x^2+1)/(3*nc^2*theta+x^2);
    f = a*K^2./((n-1).*(n-2)+K^2)+s*(n-1).*(n-2)./((n-1).*(n-2)+K^2);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    nS(i) = mean(n(j));
    
    % Hill
    H = 3;
    K = nc*((H+1)/(H-1))^(1/H);
    s = 16*H*nc/(H^2-1)/((H^2-1)*theta+4);
    a = nc*(H-1)*((H+1)^2*(theta+h)+4)/(H+1)/((H^2-1)*theta+4);
    f = a+s*n.^H./(n.^H+K^H);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    j = find(diff(sign(diff(p)))<0)+1; % maxima
    nH(i) = mean(n(j));
end

% order parameter
mS = (nS-nc)/nc;
mH = (nH-nc)/nc;

% plot inset
a1 = get(fh,'currentaxes');
a2 = axes('pos',[.575 .132 .13*[.8 1]]);
hCi = plot(h2,mS,'b-',h2,mH,'r--');
set(hCi,'linewidth',lw)
xlim([min(h2) max(h2)])
ylim([-1 1])
xlabel('$h$','fontsize',fs,'interpreter','latex')
ylabel('$m$','fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2,'xaxislocation','top',...
    'xtick',[-.2 0 .2])


print(gcf,'-depsc',[fd 'fig6.eps'])

