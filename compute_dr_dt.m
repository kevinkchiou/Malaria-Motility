function [dr2 dr dx dy] = compute_dr_dt(datastruct,ti,tf,fig)
%COMPUTE_DR_DT computes dr^2 for each track between time points ti and tf
%for a single experiment (ie datastruct(i)). It takes data structures from
%import_csv_tracks() subroutine and aggregates all the experiments in
%datastruct. If fig variable is specified, it will also plot the
%probability distribution function.
%   Detailed explanation: a data structure "datastruct" from subroutine
%import_csv_tracks() is the primary input. Specifying initial and final
%times ti and tf respectively will give displacements over that interval.
%individual track positions (datastruct(i).x(t),datastruct(i).y(t)) are
%aggregated. If the last argument is specified, a probability distribution
%of the rms displacement will be plotted with 20 bins (fixed for now).

ti_idx=zeros(numel(datastruct),1);tf_idx=zeros(numel(datastruct),1);
for i=1:numel(datastruct)
    [~,ti_idx(i)]=min(abs(datastruct(i).t-ti));
    [~,tf_idx(i)]=min(abs(datastruct(i).t-tf));
end

%vector of particle positions at r_i and r_f
xi_vec=[];yi_vec=[];xf_vec=[];yf_vec=[];
for i=1:numel(datastruct)
    dat=datastruct(i);
    xi_vec=[xi_vec,dat.x(ti_idx(i),:)];yi_vec=[yi_vec,dat.y(ti_idx(i),:)];
    xf_vec=[xf_vec,dat.x(tf_idx(i),:)];yf_vec=[yf_vec,dat.y(tf_idx(i),:)];
end

nanix=(isnan(xf_vec) | isnan(xi_vec));naniy=(isnan(yf_vec) | isnan(yi_vec));
if sum(nanix==naniy)~=numel(nanix)
	printf('Error -- compute_dr_dt(): Spatial coordinate NaN checksum error!\n');
	printf('Error -- compute_dr_dt(): Spatial coordinates disagree. Check data import routines and raw data files for inconsistencies!\n');
    error('Error -- compute_dr_dt(): Aborting...\n');
else
    nanindex=nanix;
end
rbm=~nanindex; %real bitmask, ie exclude NaN values
dx=xf_vec(rbm) - xi_vec(rbm);dy=yf_vec(rbm) - yi_vec(rbm);
dr2 = dx.^2+dy.^2;
dr = sqrt(dr2);

if nargin>3
    figure(fig);close(fig);figure(fig);
    hist(dr,20);
end
end
