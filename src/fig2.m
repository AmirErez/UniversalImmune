clear all

fh = figure(1); clf
dd2 = '../dat/';
fd = './';
fs = 16; fs2 = 12; fs3 = 14; fs4 = 20;
lw = 1.5; lw2 = 1;
ms = 6;
dx = .44; dy = .475;

nc = 500;
N = 2.5*nc;
n = (0:N)';
Z = 100;


% ---
% A: alpha (C vs theta at h = 0)
% ---

h = 0;
theta2 = linspace(-.2,.2,Z);

for i = 1:Z
    theta = theta2(i);
    
    % Schlogl: entropy and heat capacity
    x = 2*nc-3;
    K = sqrt(3*x^2+1)/2;
    s = (3*nc^3*(theta+h)+nc*x^2+x^3)/(3*nc^2*theta+x^2);
    a = ((3*x^2+1)*(3*nc^3*(theta+h)+nc*x^2+x^3)-4*x^5)...
        /(3*x^2+1)/(3*nc^2*theta+x^2);
    f = a*K^2./((n-1).*(n-2)+K^2)+s*(n-1).*(n-2)./((n-1).*(n-2)+K^2);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    nz = find(p > 0);
    SS(i) = -sum(p(nz).*log(p(nz)));
    f3 = -6*(s-a)*K^2*(3/(4*K^2-1))^(5/2);
    phi = f3*nc^2/2*[0; cumsum((n(2:end)-nc)./f(2:end))];
    CS(i) = -(1+theta)*sum(p.*(1+log(p)).*(phi - p'*phi));
    
    % Hill: entropy and heat capacity
    H = 3;
    K = nc*((H+1)/(H-1))^(1/H);
    s = 16*H*nc/(H^2-1)/((H^2-1)*theta+4);
    a = nc*(H-1)*((H+1)^2*(theta+h)+4)/(H+1)/((H^2-1)*theta+4);
    f = a+s*n.^H./(n.^H+K^H);
    c = [1; cumprod(f(2:end)./n(2:end))];
    p = c/sum(c);
    nz = find(p > 0);
    SH(i) = -sum(p(nz).*log(p(nz)));
    f3 = -(H^2-1)^2*s/8/H/nc^3;
    phi = f3*nc^2/2*[0; cumsum((n(2:end)-nc)./f(2:end))];
    CH(i) = -(1+theta)*sum(p.*(1+log(p)).*(phi - p'*phi));
end
minCS = min(CS);
minCH = min(CH);

% plot
dtheta = .03;
subplot(2,2,1)
hD = plot([-dtheta dtheta],[minCS minCS],'k-',...
    [-dtheta dtheta],[minCH minCH],'k-',...
    theta2,CS,'b-',theta2,CH,'r-');
set(hD,'linewidth',lw)
xlim([-.1 .2])
ylim([-11 5])
xlabel('Reduced temperature, $\theta$',...
    'fontsize',fs,'interpreter','latex')
ylabel('Heat capacity, $C/k_{\rm B}$',...
    'fontsize',fs,'interpreter','latex')
legend(hD(3:4),{'Schl\"ogl','Hill'},...
    'fontsize',fs2,'interpreter','latex',...
    'location','ne')
text(.005,-7,'$\alpha = 0$','horizontalalignment','center',...
    'fontsize',fs2,'interpreter','latex')
set(gca,'fontsize',fs2)
text(.14,-9,'(a)','fontsize',fs4)

% ---
% B: theta* vs 1/nc
% ---

% load data
load([dd2 'fig2_data1.mat'])

% plot
subplot(2,2,2)
hDi = plot(1./ncs,thetastarS,'b-',1./ncs,thetastarH,'r-');
set(hDi,'linewidth',lw)
xlim([0 .01])
ylim([-.01 0])
xlabel('Inverse system size, $1/n_c$','fontsize',fs,'interpreter','latex')
ylabel('Location of min., $\theta^*$',...
    'fontsize',fs,'interpreter','latex')
legend(hDi,{'Schl\"ogl','Hill'},...
    'fontsize',fs2,'interpreter','latex',...
    'location','e')
set(gca,'fontsize',fs2)
text(6e-4,-8.75e-3,'(b)','fontsize',fs4)


print(gcf,'-depsc',[fd 'fig2.eps'])

