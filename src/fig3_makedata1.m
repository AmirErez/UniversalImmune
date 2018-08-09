clear all

dd = '../dat/';
dd2 = '../dat/';

% free parameters
W = 25; % SG window (must be odd)
L = 100; % number of l bins

% fixed parameter (minimum needed for derivatives)
J = 4; % SG polynomial degree

% cosmetic parameter (for plotting only)
NI = 200; % number of I bins

% ---
% analyze data
% ---

% loop over one MEK and one SRC experiment, and drugs
expts = {'OT1Sig_20140919'};
w = 1;
for u = 1:length(expts)
    expt = expts{u}
    load([dd expt '.mat'])
    drugs = expdata.Drugs;
    for v = 1:length(drugs)
        drug = drugs{v};
        
        % save name and number of doses
        names{w} = [expt(8:end) ' ' drug];
        D(w,1) = size(expdata.RawData,2);
        
        % loop over doses
        for i = 1:D(w)
            
            % load I
            Ivals = expdata.RawData{v,i};
            S = length(Ivals); % number of samples
            
            % PDF of I
            I = linspace(min(Ivals),max(Ivals),NI)';
            dI = I(2)-I(1);
            ct = hist(Ivals,I)';
            pI = ct/S/dI;
            
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
            
            % l_c: max of phi'-phi
            phi10 = phi1 - phi0;
            j = find(~isinf(phi10));
            [ig,k] = max(phi10(j));
            jc = j(k);
            lc = l0(jc);

            % I_c, theta, h
            Ic = exp(lc);
            theta = -2*phi1(jc)/(phi1(jc)-phi3(jc));
            h = 2*phi0(jc)/(phi1(jc)-phi3(jc));
            
            % record
            Is{w}(:,i) = I;
            pIs{w}(:,i) = pI;
            ls{w}(:,i) = l;
            qs{w}(:,i) = q;
            l0s{w}(:,i) = l0;
            phi0s{w}(:,i) = phi0;
            Ics{w}(i) = Ic;
            thetas{w}(i) = theta;
            hs{w}(i) = h;
            
        end
        w = w+1;
    end
end

save([dd2 'fig3_data1.mat'])

