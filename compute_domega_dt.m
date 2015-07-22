function [datastruct_new] = compute_domega_dt(datastruct,runtracks)
%COMPUTE_DOMEGA_DT calculates angular velocity of tracks from a track data
%structure over a small time-averaged window.
%   Detailed explanation: take in a track structure output from the import
%   function import_csv_tracks, while runtracks is the output from
%   contiguous_track_stats that locates contiguous data to average over

datastruct_new=repmat(struct('t',[],'x',[],'y',[],'v',[],'r',[],'w',[],'name',''),numel(datastruct),1);
tottracks=0;
for i=1:numel(datastruct)
    tottracks=tottracks+numel(datastruct(i).x(1,:));
end
thresh = 1e5; %temporary threshold radius - roughly size of the image
for i=1:numel(datastruct)
    dat=datastruct(i);
    runs=runtracks(i).ctrack;
    len=length(dat.x(1,:));
    if len~=numel(runs)
        error('Error: tracks sizes disagree!\n');
    end
    wintime=10; %number of seconds for time-averaging window
    w=NaN(size(dat.x));vel=NaN(size(dat.x));r=NaN(size(dat.x));
    %now find contiguous tracks
    for j=1:len
        temp=runs{j};
        if ~isempty(temp) %empty case taken care of 5 lines above - NaN
            b=temp(:,1);
            e=temp(:,2);
            if numel(b)>1 || numel(e)>1
                fprintf('Warning: %s has track %d! with more than one contiguous run!\n',dat.name,j);
                for kk=1:numel(b)
                    fprintf('Warning: old (start,finish) - (%d,%d)...\n',b(kk),e(kk));
                end
                %(partially) handle contiguous run merging here
                [b e]=concatenate_contiguous_runs(b,e);
                for kk=1:numel(b)
                    fprintf('Warning: circle-fitted (start,finish) - (%d,%d)!\n',b(kk),e(kk));
                end
            end
            t=dat.t(b:e);x=dat.x(b:e,j);y=dat.y(b:e,j);
            %computes the expected range of indices given a window time wintime
            dn=floor(wintime / (mean(diff(t))));
            lidx=1+dn;ridx=e-b+1-dn;
            for k=lidx:ridx
                %find start and stops
                strt=k-dn;stp=k+dn;
                %use only subset of coordinates for circular fit
                fitx=x(strt:stp);fity=y(strt:stp);
                %use only relevant coords - exclude any NaNs that appear
                relx=fitx(~isnan(fitx));rely=fity(~isnan(fity));
                %circle_fit is sometimes "rank deficient." This happens
                %when there are tracks that are "paused"
                if mean(diff(relx))>0.01 || mean(diff(rely))>0.01
                    circ_paras=circle_fit([relx,rely]);
                    radius=circ_paras(1);%x0=circ_paras(2);y0=circ_paras(3);
                else
                    radius=0.01;
                end
                r(b+k-1,j)=radius;
                %if fit radius is too large
                if radius<thresh && radius>1
                    ttmp=t(k-1:k+1);xtmp=x(k-1:k+1);ytmp=y(k-1:k+1);
                    %since we've remapped b -> 1, the indexing requires a -1
                    [w(b+k-1,j) vel(b+k-1,j)]=compute_angular_vel(ttmp,xtmp,ytmp,circ_paras);
                else %still compute the linear velocity
                    ttmp=t(k-1:k+1);xtmp=x(k-1:k+1);ytmp=y(k-1:k+1);
                    w(b+k-1,j)=0;
                    vel(b+k-1,j)=compute_vel(ttmp,xtmp,ytmp);
                end
            end
        end
    end
    datastruct_new(i).t=dat.t;datastruct_new(i).x=dat.x;datastruct_new(i).y=dat.y;
    datastruct_new(i).w=w;datastruct_new(i).name=dat.name;
    datastruct_new(i).v=vel;datastruct_new(i).r=r;
    %these are all the tracks between the initial and final time
end
end

function [omega dvel]=compute_angular_vel(t,x,y,circ_paras)
%t,x,y should be vectors with 3 elements - middle element should be the
%evaluated point while others are to extract the velocity
r=circ_paras(1);x0=circ_paras(2);y0=circ_paras(3);
dv=[x(end)-x(1),y(end)-y(1)]/(t(end)-t(1));mid=ceil((numel(t)-1)/2)+1;
dr=[x(mid),y(mid)]-[x0,y0];drhat=dr/norm(dr);
dvperp=dv-drhat*(drhat*dv');
omega=norm(dvperp)/r;
dvel=norm(dv);
end

function dv=compute_vel(t,x,y)
dv=norm([x(end)-x(1),y(end)-y(1)]/(t(end)-t(1)));
end

function [bvec evec]=concatenate_contiguous_runs(beginvec,endvec)
%this function solves contiguous run starting from the beginning of the
%track and will truncate the rest even if there still exists more than one
%long run.
if numel(beginvec)~=numel(endvec)
    error('Error: impossible size discrepancy!\n');
end

b=beginvec(1);e=endvec(1);

%first two are stopping conditions, third is recursion.
if length(beginvec)==1
    bvec=beginvec;evec=endvec;
elseif beginvec(2) - endvec(1) >2
    bvec=beginvec;evec=endvec;
elseif beginvec(2) - endvec(1) < 3
    beginvec=[beginvec(1);beginvec(3:end)];
    endvec=endvec(2:end);
    [bvec evec]=concatenate_contiguous_runs(beginvec,endvec);
else
    error('Error: recursion failure.\n');
end

end
