function [x xinv pr prinv rinv_agg,r_agg] = plot_circ_probdistfcn(cell_datastruct,nbin,figs)
%PLOT_CIRC_PROBDISTFCN Takes outputs from a csv import script and creates the probability density / distribution of curvatures. It requires the extra scripting functions CONTIGUOUS_TRACK_STATS, COMPUTE_DOMEGA_DT, and AGGREGATE_CIRC_DATA. The arguments are a cell of datastructures from all tracks, the number of bins, and two figure handles: one for the radius and the other for the curvature. It should be noted that CONTIGUOUS_TRACK_STATS also requires the CONTIGUOUS and the CIRCLE_FIT subroutines
%   This function plots circle and curvature probability distributions of an aggregation of tracked sporozoites. There are fail-safes written in to handle non-contiguous tracks by splitting them up, or joining them together in case of short breaks. This requires five other custom-written subroutines to function, as detailed in the description.
if nargin<2
    nbin=40;
end
if nargin<3
    figs=[1 2];
end
num=numel(cell_datastruct);
leg=cell(num,1);
rinv=cell(num,1);r=cell(num,1);
nrinv=zeros(num,nbin);nr=zeros(num,nbin);prinv=zeros(num,nbin);pr=zeros(num,nbin);
rinv_agg=[];r_agg=[];
for i=1:num
    dat=cell_datastruct{i};
    leg{i}=dat(1).name(1:end-7);
    if ~isfield(dat,'r') %somehow there's no r-field
        %find the contiguous track runs
        [~,temp_runtracks]=contiguous_track_stats(dat);
        %compute things like radius, etc
        dat=compute_domega_dt(dat,temp_runtracks);
    end
    [rinv{i} r{i}]=aggregate_circ_data(dat);
    rinv_agg=[rinv_agg;rinv{i}];r_agg=[r_agg;r{i}];
end

%eliminates inferred radii that are larger than the physical window size.
r_agg=r_agg(r_agg<1e4);

[~,xinv]=hist(rinv_agg,nbin);
[~,x]=hist(r_agg,nbin);
for i=1:num
    nrinv(i,:)=hist(rinv{i},xinv);
    prinv(i,:)=nrinv(i,:)/sum(nrinv(i,:));
    nr(i,:)=hist(r{i},x);
    pr(i,:)=nr(i,:)/sum(nr(i,:));
end

figure(figs(1));plot(xinv,prinv);xlabel('Curvature (\kappa)');ylabel('P(\kappa)');legend(leg,'Location','NorthEast');
figure(figs(2));plot(x,pr);xlabel('Radius (r)');ylabel('P(r)');legend(leg,'Location','NorthEast');
end

