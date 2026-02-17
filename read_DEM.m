
function [minlat maxlat minlon maxlon reslat reslon elev] = read_DEM(fileName)

[filepath,name,ext] = fileparts(fileName);
if ext == ".bin"

    % Read 16 bit SRTM30+ file, and keep it int16 to save memory.
    fid = fopen(fileName,'r');
    %[topography, cnt] = fread(fid, inf, 'int16=>int16');
    [data, cnt] = fread(fid, inf, 'double','native');
    
    fclose(fid);
    minlat=data(1);
    maxlat=data(2);
    minlon=data(3);
    maxlon=data(4);
    reslat=data(5);
    reslon=data(6);
    % make image a rectangle. SRTM30+ has 4800 columns north of 60S,
    % but 7200 columns south of 60S. So user has say how many. . .
    numCols=(maxlon-minlon)/reslon+1;
    numRows = (cnt-6) / numCols;
    elev=data(7:end);
    elev = reshape(elev,round(numCols),round(numRows))';
else
    [A,R] = readgeoraster(fileName);
    minlat = R.LatitudeLimits(1);
    maxlat = R.LatitudeLimits(2);
    minlon = R.LongitudeLimits(1);
    maxlon = R.LongitudeLimits(2);
    if class(R)== "map.rasterref.GeographicCellsReference"
        reslat = R.CellExtentInLatitude;
        reslon = R.CellExtentInLongitude;
        minlat = minlat + reslat/2;
        maxlat = maxlat - reslat/2;
        minlon = minlon + reslon/2;
        maxlon = maxlon - reslon/2;
     else
        reslat = R.SampleSpacingInLatitude;
        reslon = R.SampleSpacingInLongitude;
     end
        elev = flip(double(A));

end

end