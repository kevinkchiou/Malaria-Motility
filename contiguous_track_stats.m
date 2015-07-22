function [runlengths runcells datastruct_out] = contiguous_track_stats(datastruct)
%CONTIGUOUS_TRACK_STATS Takes a data structure output by IMPORT_CSV_TRACKS and finds contiguous sets of data within each track. If the breaks are short and below threshold, it will join the tracks. Otherwise it will separate them out into separate runs and create a new datastructure in the form of "datastruct_out"
%   This function requires the subroutine CONTIGUOUS which was written for more general purposes. Much of the complexity here is dealing with the specific datastructures output by IMPORT_CSV_TRACKS

runlengths = cell(numel(datastruct),1);
runcells=repmat(struct('ctrack',[]),numel(datastruct),1);
datastruct_out=datastruct;
for i=1:numel(datastruct)
    dat=datastruct(i);
    if(sum(sum((isnan(dat.x)==isnan(dat.y))))~=numel(dat.x))
        %check for zeroth order inconsistency in tracks
		%often happens due to comma-type or two/three
		%column data differences
        error('Error: inconsistent coordinates');
    else
		%use a simple binary type of matrix for simpler
		%contiguous computation
        mat=~isnan(dat.x);
        matconj=~mat;
    end
    
    %use contiguous.m function written separately
    sz=size(mat);
    %some corrective measures for small breaks in positional information
    for j=1:sz(2)
        %check that there are terms with NaN
        if sum(matconj(:,j)==1)>0
            rtemp=contiguous(matconj(:,j),1);
            runlen=rtemp{1,2}(:,2)-rtemp{1,2}(:,1)+1;
            %if NaN sequence is short, then interpolate the data
            shortrunidx=runlen<3;
            previdx=rtemp{1,2}(shortrunidx,1)-1;nextidx=rtemp{1,2}(shortrunidx,2)+1;
            for k=1:numel(previdx)
                if previdx(k)<1
                    %flatten the profile
                    dat.x(previdx(k)+1:nextidx(k)-1,j) = dat.x(nextidx(k),j);
                    dat.y(previdx(k)+1:nextidx(k)-1,j) = dat.y(nextidx(k),j);
                elseif nextidx>numel(matconj(:,j))
                    %flatten the profile
                    dat.x(previdx(k)+1:nextidx(k)-1,j) = dat.x(previdx(k),j);
                    dat.y(previdx(k)+1:nextidx(k)-1,j) = dat.y(previdx(k),j);
                else
                    %interpolate (linear) the profile
                    tvec=dat.t(previdx(k):nextidx(k));dt=tvec(end)-tvec(1);
                    vx=dat.x(nextidx(k),j)-dat.x(previdx(k),j)/dt;
                    vy=dat.y(nextidx(k),j)-dat.y(previdx(k),j)/dt;
                    dx=dat.x(previdx(k),j)+vx*diff(tvec);
                    dy=dat.y(previdx(k),j)+vy*diff(tvec);
                    dat.x(previdx(k)+1:nextidx(k)-1,j) = dx(1:end-1);
                    dat.y(previdx(k)+1:nextidx(k)-1,j) = dy(1:end-1);
                end
            end
        end
    end
    datastruct_out(i)=dat;
    %Aggregate the run lengths for each dataset
    temprunlen=[];ctrack_temp=cell(sz(2),1);
    for j=1:sz(2)
        %check to make sure there are real elements
        if sum(mat(:,j)==1)>0
            %since I am passing a scalar into the second argument of
            %contiguous(), the resulting cell is of size {1,2} for the one
            %contiguous number, and its start and stop matrix
            rtemp=contiguous(mat(:,j),1);
            temprunlen=[temprunlen;rtemp{1,2}(:,2)-rtemp{1,2}(:,1)+1];
            ctrack_temp{j}=[rtemp{1,2}(:,1),rtemp{1,2}(:,2)];
        end
    end
    runlengths{i}=temprunlen;
    runcells(i).ctrack=ctrack_temp;
end
end
