function [allSSEs, I1s, npoints, fitresult] = fitCv(delta_theta, D)
dd3 = '../Cvdat/';
files = glob([dd3 'figCv_data1_I*.mat']);
I1s = nan(length(files),1);
SSEs = zeros(length(files),length(D));

for r = 1:length(I1s)
    SSE_cnt = 0;
%     disp(r)
%     load([dd3 'figCv_data1_I' num2str(r) '.mat'])
    loaded = load(files{r});
    I1s(r) = loaded.I1;
    disp(I1s(r));
    for w = 1:length(loaded.D)
        SSE = 0;
        cnt = 0;
        for i = 1:D(w)
            if(abs(loaded.thetas{w}(i))>delta_theta), continue; end
            theory = interp1(loaded.tH, loaded.CH, loaded.thetas{w}(i));
            if(~isnan(theory))
                SSE = SSE + (theory - loaded.Cs{w}(i))^2;
                cnt = cnt + 1;
                SSE_cnt = SSE_cnt +1;
            end
%             plot(thetas{w}(i),Cs{w}(i),sh(w),'markersize',ms(i),...
%                 'markerfacecolor',co(j,:),'markeredgecolor','k');
        end
        SSEs(r,w) = SSE/cnt;
    end
    disp(['I1:' num2str(I1s(r)) ' in SSE:' num2str(SSE_cnt)]);
end

disp(['Delta theta ' num2str(delta_theta) ' ; using ' num2str(SSE_cnt) ' fitting points']);
npoints = SSE_cnt;
allSSEs = nansum(SSEs,2);
ft = fittype( 'A + ((x-x0)^2)/(2*sigma^2)', 'independent', 'x', 'dependent', 'y' );
% exclude_inds = [find(I1s<0.2),find(I1s>0.8)];
% excludedPoints = excludedata( I1s', allSSEs, 'Indices', exclude_inds );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [min(allSSEs) 0.1 0.1];
% opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( I1s, allSSEs, ft, opts );
