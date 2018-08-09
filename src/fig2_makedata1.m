clear all
dd2 = '../dat/';

Z1 = 20;
Z2 = 2000;
nc = 500;
N = 2.5*nc;
n = (0:N)';
h = 0;

% data for Fig 1d inset: theta^* vs 1/n_c
ncs = 1./linspace(.001,.01,Z1);
theta3 = linspace(-.1,0,Z2);
theta3_ = theta3(1:end-1);
for j = 1:Z1
    nc = ncs(j)
    N = 3*nc;
    n = (0:N)';
    SSi = []; SHi = []; CSi = []; CHi = [];
    for i = 1:Z2
        theta = theta3(i);
        
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
        SSi(i) = -sum(p(nz).*log(p(nz)));
        f3 = -6*(s-a)*K^2*(3/(4*K^2-1))^(5/2);
        phi = f3*nc^2/2*[0; cumsum((n(2:end)-nc)./f(2:end))];
        CSi(i) = -(1+theta)*sum(p.*(1+log(p)).*(phi - p'*phi));
        
        % Hill: entropy and heat capacity
        H = 3;
        K = nc*((H+1)/(H-1))^(1/H);
        s = 16*H*nc/(H^2-1)/((H^2-1)*theta+4);
        a = nc*(H-1)*((H+1)^2*(theta+h)+4)/(H+1)/((H^2-1)*theta+4);
        f = a+s*n.^H./(n.^H+K^H);
        c = [1; cumprod(f(2:end)./n(2:end))];
        p = c/sum(c);
        nz = find(p > 0);
        SHi(i) = -sum(p(nz).*log(p(nz)));
        f3 = -(H^2-1)^2*s/8/H/nc^3;
        phi = f3*nc^2/2*[0; cumsum((n(2:end)-nc)./f(2:end))];
        CHi(i) = -(1+theta)*sum(p.*(1+log(p)).*(phi - p'*phi));
    end
    
    % theta^*
    [minCSi,iS] = min(CSi);
    [minCHi,iH] = min(CHi);
    thetastarS(j) = theta3(iS);
    thetastarH(j) = theta3(iH);
end

save([dd2 'fig2_data1.mat'])