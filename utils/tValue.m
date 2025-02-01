function [mn,se]= tValue(gmeans,gcov)
% Make sure NaN groups don't affect other results
t = isnan(gmeans);
if any(t)
    gcov(t,:) = 0;
    gcov(:,t) = 0;
end
ng = length(gmeans);
M = nchoosek(1:ng, 2);      % all pairs of group numbers
g1 = M(:,1);
g2 = M(:,2);
mn = gmeans(g1) - gmeans(g2);
i12 = sub2ind(size(gcov), g1, g2);
gvar = diag(gcov);
se = sqrt(gvar(g1) + gvar(g2) - 2 * gcov(i12));              

