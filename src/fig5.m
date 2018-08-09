clear all

dd = '../dat/fig5/';
fd = './';
fs = 16; fs2 = 12; fs3 = 20;
lw = 1;
ms = 5;
gr = [0 .5 0];

param = load([dd 'multipass3.param.dat']);
thetamin = param(1);
thetamax = param(2);
h = param(3);
nc = param(4);
lrhomin = param(5);
lrhomax = param(6);
T = param(7);
Ntheta = param(8);
Nrho = param(9);
Nx = param(10);

x = (0:Nx)';
rhos = logspace(lrhomin,lrhomax,Nrho);
thetas = linspace(thetamin,thetamax,Ntheta)';

% SG filter for S(theta)
W1 = 5;
J1 = 1;
[B1,G1] = sgolay(J1,W1);
Y1 = (W1+1)/2;
dtheta = thetas(2)-thetas(1);
theta0 = thetas(Y1:Ntheta-Y1);

% SG filter for f(x)
W2 = 25;
J2 = 3;
[B2,G2] = sgolay(J2,W2);
Y2 = (W2+1)/2;
x0 = x(Y2:Nx-Y2);

% for inferring nc
j1 = round(.8*nc);
j2 = round(1.2*nc);


for l = 1:Nrho
    rho = rhos(l);
    lstr{l} = ['$\rho = ' num2str(rho) '$'];
    
    for k = 1:Ntheta
        theta = thetas(k);

        c = load([dd 'multipass3_rho' num2str(l-1) ...
            '_theta' num2str(k-1) '.Cx.dat']);
        p = c/T;

        % ---
        % 1. analytic theta from deterministic eq,
        %    C using entropy and numerical derivative
        % ---
        nz = find(p>0);
        S(k,l) = -p(nz)'*log(p(nz));
        
        % ---
        % 2. inferred theta from coarse-grained model f(x),
        %    C using associated expression involving f
        % ---

        f = [0; x(2:end).*p(2:end)./p(1:end-1)];
        df = gradient(f);
        d2f = gradient(df);
        d3f = gradient(d2f);

        % SG filter
        for j = Y2:Nx-Y2
            f0(j-Y2+1) = G2(:,1)'*f(j-Y2+1:j+Y2-1);
            f1(j-Y2+1) = G2(:,2)'*f(j-Y2+1:j+Y2-1);
            f2(j-Y2+1) = G2(:,3)'*f(j-Y2+1:j+Y2-1);
            f3(j-Y2+1) = G2(:,4)'*f(j-Y2+1:j+Y2-1);
        end

        % infer Ising parameters
        [ig,j_] = max(f1(j1-Y2:j2-Y2));
        j = j_ + j1-Y2 - 1;
        nc_ = x0(j);
        fc = f0(j);
        f1c = f1(j);
        f3c = f3(j);
        theta_(k,l) = -2*(1-f1c)/f3c/nc_^2;

        % calculate C
        nz = find(p > 0);
        f(find(isnan(f))) = inf;
        phi = f3c*nc_^2/2*[0; cumsum((x(2:end)-nc_)./f(2:end))];
        C_(k,l) = -(1+theta_(k,l))*...
            sum(p(nz).*(1+log(p(nz))).*(phi(nz)-p(nz)'*phi(nz)));
    end
    
    % SG filter
    for j = Y1:Ntheta-Y1
        S0(j-Y1+1,l) = G1(:,1)'*S(j-Y1+1:j+Y1-1,l);
        S1(j-Y1+1,l) = G1(:,2)'*S(j-Y1+1:j+Y1-1,l)/dtheta;
    end
    C(:,l) = (theta0+1).*S1(:,l);
end

% plot
figure(1); clf
subplot(2,2,1)
hA = plot(theta0,C,'o-');
set(hA,'linewidth',lw,'markersize',ms)
set(hA(1),'color','b','marker','o')
set(hA(2),'color','r','marker','s')
set(hA(3),'color',gr,'marker','^')
ylim([-5 5])
xlabel('Reduced temperature, $\theta$',...
    'fontsize',fs,'interpreter','latex')
ylabel('Heat capacity, $C/k_{\rm B}$',...
    'fontsize',fs,'interpreter','latex')
legend(lstr,'location','ne',...
    'fontsize',fs2,'interpreter','latex')
% title('Exact',...
%     'fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2)
text(-.09,-4,'(a)','fontsize',fs3)

subplot(2,2,2)
hB = plot(theta_,C_,'o');
set(hB,'linewidth',lw,'markersize',ms)
set(hB(1),'color','b','marker','o')
set(hB(2),'color','r','marker','s')
set(hB(3),'color',gr,'marker','^')
xlim([-1 2.5])
ylim([-2 2])
xlabel('Reduced temperature, $\theta$',...
    'fontsize',fs,'interpreter','latex')
ylabel('Heat capacity, $C/k_{\rm B}$',...
    'fontsize',fs,'interpreter','latex')
legend(lstr,'location','ne',...
    'fontsize',fs2,'interpreter','latex')
% title('Coarse-grained',...
%     'fontsize',fs,'interpreter','latex')
set(gca,'fontsize',fs2)
text(-.87,-1.6,'(b)','fontsize',fs3)

print(gcf,'-depsc',[fd 'fig5.eps'])

