function data = import_csv_tracks(FILESTR,DIRNAME,option,commaoption)
%IMPORT_CSV_TRACKS Takes .csv files of tracks from the cell of regex strings FILESTR in directory DIRNAME and creates a data structure for each element in FILESTR. "option" is a string detailing whether the data is two or three column output from the tracker (ie (x,y) or (x,y,z)). commaoption is whether or not the lines end with a comma
%   This subroutine imports tracks from a set of .csv files. Much of the complication is dealing with various possible outputs from the tracking software. Otherwise the standard readcsv from MATLAB would work just as well.

if nargin<2
    DIRNAME='./';
end
if nargin<3
    option='TwoColumn';
end
if nargin<4
    commaoption=0;
end

len=numel(NAMES);
data=repmat(struct('t',[],'x',[],'y',[],'name',''),len,1);
for i=1:len %len-1 for last file with empty slots
    currname=NAMES(i).name;
    %fprintf('Getting data from %s...\n',currname); %for debugging
    dat=cat_csvread([DIRNAME,'/',currname],commaoption); %my more robust csvread
    %function here to extract tracks correctly (x_i,y_i,z_i) vs (x_i,y_i)
    %option = 'ThreeColumn'; %option = 'TwoColumn';
    data(i) = extract_tracks_from_matrix(dat,currname,option);
end
end

function out=cat_csvread(dirstr,commaopt)
%get number of rows in file
[Nrow,Nelt,Srow]=extract_file_info(dirstr,commaopt);
fptr=fopen(dirstr,'r');
tempbuffer=find_first_data_row(fptr);
count=0;out=zeros(Nrow-Srow,Nelt);
while ~feof(fptr)
    count=count+1;
    if tempbuffer(end)~=','
        buff=[tempbuffer,','];
    elseif commaopt~=0
        buff=[tempbuffer,','];
    else
        buff=tempbuffer;
    end
    c_idx=strfind(buff,','); %comma locations
    d1_idx=diff([0,c_idx]); %distances between commas
    d2_idx=diff(c_idx); %shifted distances between commas
    endlocs=c_idx(d1_idx>1)-1; %locs with relevant distances between commas
    startlocs=[1,c_idx(d2_idx>1)+1]; %create startloc vector from endlocs
    vals=NaN(size(c_idx));j=0;
    for i=1:numel(vals)
        if d1_idx(i)>1
            j=j+1;
            vals(i)=str2double(buff(startlocs(j):endlocs(j)));
        end
    end
    if numel(out(count,:))~=numel(vals)
        fprintf('Error: near line %d, in file %s. Values follow:\n',count,dirstr);
        vals
        size(out(count,:))
    end
    out(count,:) = vals;
    %now we have extracted the values, and we have to find the track
    %that each of them goes into
    tempbuffer=fgetl(fptr);
end
%out=zeros(row+count,length(c_idx)+1); %temporarily set zero matrix
fclose(fptr);
end

function [Nrow,Nelt,Skiprow]=extract_file_info(dirstr,buffopt)
fptr=fopen(dirstr,'r');
[buff,Skiprow]=find_first_data_row(fptr);
Nrow=Skiprow;
if buffopt~=0
    buff=[buff,','];
end
commalocs=strfind(buff,',');
if strcmp(buff(end),',')
    Nelt=numel(commalocs);
else
    Nelt=numel(commalocs)+1;
end
while ~feof(fptr)
    fgetl(fptr);Nrow=Nrow+1;
end
fclose(fptr);
end

%recursive solution to finding the first row of data
function [outbuffer,n]=find_first_data_row(fptr,n)
if nargin<2
    n=0;
else
    n=n+1;
end
outbuffer=fgetl(fptr);
if numel(outbuffer)==0 %line with no characters
    outbuffer=fgetl(fptr);n=n+1;
end
if strcmp(outbuffer(1),',')
    fprintf('Warning: the line buffer starts with a comma!\n');
end
if isempty(outbuffer)
    outbuffer=fgetl(fptr);n=n+1;
end
while isnan(str2double(outbuffer(1)))
    [outbuffer,n]=find_first_data_row(fptr,n);
end
end

%takes in a matrix of values from csv importing and associates with tracks
function out=extract_tracks_from_matrix(in,name,option)
out=struct('t',[],'x',[],'y',[],'name','');
out.t=in(:,1)/1000;
out.name=name;
%zeroth order check where if the number of spatial coordinates plus t comes
%to an even number, then it must be three column.
sz=size(in);
if strcmp(option,'TwoColumn') && mod(sz(2),2)==0
    option='ThreeColumn';
    fprintf('Warning: The data format for %s does not match two-column!\n',name);
end
if strcmp(option,'TwoColumn') %other test for possible 3-column
    coord1=in(:,2:3:end);
    coord2=in(:,3:3:end);
    coord3=in(:,4:3:end);
    dd1=var(mean(coord1,1));dd2=var(mean(coord2,1));dd3=var(mean(coord3,1));
    if (dd1<1e-3 || dd2<1e-3 || dd3<1e-3)
        fprintf('Warning: likely a 3-column format!\nWarning: Extremely ');
        fprintf('low variance coordinate in file %s:\n',name);
        fprintf('Warning: (v1,v2,v3)=(%f,%f,%f)\n',dd1,dd2,dd3);
    end
end
switch option
    case 'TwoColumn'
        out.x=in(:,2:2:end); %pick out x and y data by skipping
        out.y=in(:,3:2:end);
    case 'ThreeColumn'
        sz1=size(in(:,2:3:end));
        sz2=size(in(:,3:3:end));
        sz3=size(in(:,4:3:end));
        if sz1(2)~=sz2(2) || sz1(2)~=sz3(2) || sz2(2)~=sz3(2)
            fprintf('Error: file %s\n',name);
            error('Error: Maybe not 3-column format!');
        else
            num=sz1(2);
        end
        
        %performed this slower way instead of vectorized due to z coord
        %actually getting mixed up in order compared to x and y.
        for i=1:num
            dd1=var(in(:,2+(i-1)*3));dd2=var(in(:,3+(i-1)*3));dd3=var(in(:,4+(i-1)*3));
            [~,idx]=min([dd1,dd2,dd3]);
            idx1=mod(idx,3)+2;idx2=mod(idx+1,3)+2;
            out.x(:,i)=in(:,idx1+(i-1)*3);
            out.y(:,i)=in(:,idx2+(i-1)*3);
        end
        
        %use a slightly less efficient version for taking the mean to
        %properly maintain matrix structure;
        %m1=mean_exNaN(coord1);m2=mean_exNaN(coord2);m3=mean_exNaN(coord3);
        
        %different ways of possibly identifying relevant z-direction,
        %by comparing statistics between tracks.
        %df1=diff(m1);df2=diff(m2);df3=diff(m3);[~,idx]=min([df1,df2,df3]);
        %dd1=var(m1);dd2=var(m2);dd3=var(m3);[~,idx]=min([dd1,dd2,dd3]);
        %mm1=mean(m1);mm2=mean(m2);mm3=mean(m3);[~,idx]=min([mm1,mm2,mm3]);
        %idx1=mod(idx,3)+2;idx2=mod(idx+1,3)+2;
        %out.x=in(:,idx1:3:end);
        %out.y=in(:,idx2:3:end);
    otherwise
        error('Error: Invalid choice!\n');
end
%zeroth order check before we let you output complete nonsense.
if sum(sum(isnan(out.x)==isnan(out.y)))==numel(out.x)
    return;
else
    error('Error: zeroth order inconsistency in coordinates!\n');
end
end

function out=mean_exNaN(in,dim)
if nargin<2
    dim=1;
end
sz=size(in);
switch dim
    case 1
        in=in;out=zeros(1,sz(2));
    case 2
        in=in';out=zeros(sz(1),1);
    otherwise
        error('Cannot handle other values for dim at the moment!\n');
end
sz=size(in);
for i=1:sz(2)
    vec=in(:,i);
    out(i)=mean(vec(~isnan(vec)));
end
end
