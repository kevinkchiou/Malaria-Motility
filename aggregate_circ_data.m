function [rinv r r_raw] = aggregate_circ_data(datastruct,threshold)
%AGGREGATE_CIRC_DATA Simple way of aggregating similar curvature and radius data together while consistently eliminating spurious radius fits stemming from inference sensitivity
%   The linear approach to inferring radii has occasional sensitive results when tracks are extremely straight and somewhat noisy. This can be eliminated to some degree by using longer time windows, but can nonetheless appear. Here, during the aggregation process of similar type tracks, we choose consistent threshold sensitivity as well as eliminating tracks with their radii that were automatically set to NaN due to non-contiguous nature of some pieces of data.

if nargin<2
    threshold = 0.01;
end
r_raw=[];
for i=1:numel(datastruct)
    r_raw=[r_raw;datastruct(i).r(~isnan(datastruct(i).r))];
end
r=r_raw(r_raw>threshold);
rinv=1./r;
rinv=rinv(~isinf(rinv));
r=r(~isinf(r));
end

