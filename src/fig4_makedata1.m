clear all

dd = '../dat/';
dd2 = '../dat/';

yls = [-80 20; -80 20; -80 20; -170 20];
ys = yls(:,1) + diff(yls,1,2)/8;

% free parameters
W = 25; % SG window (must be odd)
L = 100; % number of l bins

% free parameter used only for heat capacity
I1s = [10 1 .1 .01]; % intensity of one molecule

% fixed parameter (minimum needed for derivatives)
J = 4; % SG polynomial degree

% ---
% analyze data
% ---

for r = 1:length(I1s)
    I1 = I1s(r)
    if r == 4; disp('experiment'); end

    % loop over all experiments and drugs
    expts = {'OT1Sig_20140918','OT1Sig_20140919','OT1Sig_20151206'};
    w = 1;
    for u = 1:length(expts)
        expt = expts{u};
        load([dd expt '.mat'])
        drugs = expdata.Drugs;
        for v = 1:length(drugs)
            if r == 4; disp([num2str(round(100*w/9)) '%']); end
            drug = drugs{v};
            
            % save name and number of doses
            names{w} = [expt(8:end) ' ' drug];
            D(w,1) = size(expdata.RawData,2);
            
            % loop over doses
            for i = 1:D(w)
                
                % load I data
                Ivals = expdata.RawData{v,i};
                S = length(Ivals); % number of samples
                
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
                q0 = []; phi0 = []; phi1 = []; phi2 = []; phi3 = [];
                for j = Y:L-Y
                    q0(j-Y+1,1) = G(:,1)'*q(j-Y+1:j+Y-1);
                    phi0(j-Y+1,1) = G(:,2)'*Q(j-Y+1:j+Y-1)/dl;
                    phi1(j-Y+1,1) = 2*G(:,3)'*Q(j-Y+1:j+Y-1)/dl^2;
                    phi2(j-Y+1,1) = 6*G(:,4)'*Q(j-Y+1:j+Y-1)/dl^3;
                    phi3(j-Y+1,1) = 24*G(:,5)'*Q(j-Y+1:j+Y-1)/dl^4;
                end
                l0 = l(Y:L-Y);
                
                % l_c: max of phi'-phi
                phi10 = phi1 - phi0;
                j = find(~isinf(phi10));
                [ig,k] = max(phi10(j));
                jc = j(k);
                lc = l0(jc);
                nc = exp(lc)/I1;
                
                % theta, h
                theta = -2*phi1(jc)/(phi1(jc)-phi3(jc));
                h = 2*phi0(jc)/(phi1(jc)-phi3(jc));
                
                % heat capacity
                l_ = log(I1*n(2:end));
                f = interp1(l0,phi0,l_) + n(2:end);
                f(find(isnan(f))) = inf;
                psi = (phi3(jc)-phi1(jc))/2 ...
                    *[0; cumsum((n(2:end)/nc-1)./f)];
                k = find(pn > 0);
                C = -(1+theta)*sum(pn(k).*(1+log(pn(k))) ...
                    .*(psi(k)-pn'*psi));
                
                % record
                ncs{w}(i) = nc;
                hs{w}(i) = h;
                thetas{w}(i) = theta;
                Cs{w}(i) = C;
                
            end
            w = w+1;
        end
    end
    
    if r == 4; disp('theory'); end
    xl = [-1 1.2];
    yl = yls(r,:);
    
    % find mean nc from data; use it for theory
    all_nc = []; all_C = []; all_h = [];
    for w = 1:length(D)
        for i = 1:D(w)
            all_nc = [all_nc ncs{w}(i)];
            all_C = [all_C Cs{w}(i)];
            all_h = [all_h hs{w}(i)];
        end
    end
    nc = mean(all_nc);
    
    N = 2.5*nc;
    n = (0:N)';
    H = 4;
    h = 0;
    Z = 1000;
    tmin1 = -4./(H+1).^2*.625; % where a < 0
    tH = linspace(tmin1,xl(2),Z);
    
    for i = 1:Z
        if r == 4 && mod(i,Z/10) == 0
            disp([num2str(round(100*i/Z)) '%'])
        end

        % Hill (compute p_n starting from middle, for numerical stability)
        t = tH(i);
        K = nc*((H+1)/(H-1))^(1/H);
        s = 16*H*nc/(H^2-1)/((H^2-1)*t+4);
        a = nc*(H-1)*((H+1)^2*(t+h)+4)/(H+1)/((H^2-1)*t+4);
        f = a+s*n.^H./(n.^H+K^H);
        [ig,j] = min(abs(n-nc));
        cp = cumprod(f(j+1:end)./n(j+1:end));
        cm = cumprod(n(j:-1:2)./f(j:-1:2));
        c = [cm(end:-1:1); 1; cp];
        p = c/sum(c);
        f3 = -(H^2-1)^2*s/8/H/nc^3;
        psi = f3*nc^2/2*[0; cumsum((n(2:end)-nc)./f(2:end))];
        nz = find(p > 0);
        CH(i) = -(1+t)*sum(p(nz).*(1+log(p(nz)))...
            .*(psi(nz) - p(nz)'*psi(nz)));
    end
    
    save([dd2 'fig4_data1_I' num2str(r) '.mat'])

end
