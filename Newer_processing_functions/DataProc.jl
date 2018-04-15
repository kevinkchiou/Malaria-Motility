function extractXLSFilenames()
    files = readdir("./");
    bm = [ifelse(search(x,".xls")!=0:-1,true,false) for x in files];
    return files[bm];
end

function convertDataFrameToFloat(df)
    #convert time to seconds
    import DataFrames
    df[1]=df[1]/1000.0;
    for i in 2:length(df)
        df[i]=[ifelse(z<-50,NA,z) for z in [parse(Float64,y) for y in [ifelse(DataFrames.isna(x),"-100.0",x) for x in df[i]]]]
    end
    return df
end

function timeSliceDataFrame(df,interval=3.0)
    #determine sizes of new dataframes
    times = df[1,1]:interval:df[end,1]
    numtimes = length(times)
    numtracks=Int((length(df)-1) / 3)
    positions = zeros(numtimes,2*numtracks) #convert from old times to new times
    #populate sparse posconvmat matrix
    for j in 1:numtracks,(i,t) in enumerate(times)
        aft=find(df[1]-t .>= 0)[1]
        bef=find(df[1]-t .<= 0)[end]
        if bef!=aft
            dt2 = df[aft,1] - t
            dt1 = t - df[bef,1]
            dt=dt1+dt2
            #linear interpolation for x and y in df -- exclude z
            poaitions[i,2*j-1] = (df[bef,3*j-2]*dt2 + df[aft,3*j-2]*dt1) / dt
            positions[i,2*j] = (df[bef,3*j-1]*dt2 + df[aft,3*j-1]*dt1) / dt
        else
            positions[i,2*j-1] = df[bef,3*j-2]
            positions[i,2*j] = df[bef,3*j-1]
        end
    end
    #now all times and positions of new interpolated data are created
end
