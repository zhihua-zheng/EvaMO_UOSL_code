function S = pinBin(bi,p,q,eM)
%
% pinBin
%==========================================================================
%
% USAGE:
%  S = pinBin(bi,p,q,eM)
%
% DESCRIPTION:
%  
%
% INPUT:
%
%  zi - 1-D column vector of designed bin edge postions in ascending order
%  p  - 1-D column vector of coordinates for quantity be be parsed
%  q  - 1-D column vector of quantity at according location in p
%  eM - 'left' or 'none', option for including bin edges
%
% OUTPUT:
%
%  S - struct contains statistical information for each bin
%      qm   - mean value
%      n    - number of available data
%      qstd - standard deviation
%      dof  - degree of freedom
%      tupp - upper bound for certain significance level in student's
%             distribution
%
% AUTHOR:
%  June 26 2019, Zhihua Zheng                             [ zhihua@uw.edu ]
%==========================================================================

%%

% edges of bins
% 2-column matrix containing edges of bins in ascending order. 1st column 
% specifies lower edge, while 2nd column specifies upper edge.
bins = [bi(1:end-1) bi(2:end)];
nbin = size(bins,1);
nq   = size(q,2);

% confidence level bounds
conf = 0.95;
alpL = (1-conf)/2;
alpU = 1-alpL;

% loop through each bins
S.qm   = nan(nbin,nq);      % sample mean
S.n    = zeros(size(S.qm)); % sample size
S.qstd = zeros(size(S.qm)); % standard deviation
S.dof  = zeros(size(S.qm));
S.tupp = zeros(size(S.qm)); 
S.qmL  = zeros(size(S.qm)); 
S.qmU  = zeros(size(S.qm)); 

for i = 1:nbin
    
    switch eM
        case 'left'
            bInx = p >= bins(i,1) & p < bins(i,2);
            
        case 'none'
            bInx = p > bins(i,1) & p < bins(i,2);
    end
             
    qBin        = q(bInx,:);
    S.n(i,:)    = sum(~isnan(qBin));
    S.qm(i,:)   = nanmean(qBin);
    S.qstd(i,:) = nanstd(qBin);
    
    if min(S.n(i,:)) >= 10
        
    % empirical bootstrap method for confidence limits of sample mean
    % resample size is the same as original sample
    % resampling with replacement from the rows of original data
    nboot     = 1000; % numebr of resamples
    [bootM,~] = bootstrp(nboot,@nanmean,qBin');
    
    % the variation of sample mean S.qm is well approximated by the
    % variation of resmaple mean bootm
    boot_delta  = bootM - repmat(S.qm(i,:),nboot,1);
    qm_delta    = prctile(boot_delta,100*[alpL,alpU]);
    S.qmL(i,:)  = -qm_delta(1);
    S.qmU(i,:)  =  qm_delta(2);
    
    else
    
    S.qmL(i,:) = NaN;
    S.qmU(i,:) = NaN;
        
    end
end

S.dof = S.n - 1; % degree of freedom
S.dof(S.dof <= 0) = NaN;

% critical vlaues for Student's distribution
S.tupp = tinv(alpU,S.dof);
% S.tlow = tinv(alpL,S.dof);

S.qerr = S.tupp .* S.qstd ./ sqrt(S.dof); % symmetric confidence limits

end

