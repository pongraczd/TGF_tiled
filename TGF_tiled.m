function TGF_tiled(ComPoints, tiles, DetailedDEM, DetailedREF, MassDensity, TessDEM, TessREF,TessDensity, CoarseDEM, CoarseREF,CoarseDensity,GlobalDEM,GlobalREF, outname,ikind, itype,idensity, flag_earth,rzones, e)
    output_report=strcat(outname,'_report.txt');  
    comPoints = load(ComPoints);
    lat = comPoints(:,2);
    lon = comPoints(:,3);
    num=size(comPoints,1);
    minlat = min(lat);
    maxlat = max(lat);
    minlon = min(lon);
    maxlon = max(lon);
    [DetailedDEM_raster , DetailedDEM_rasterref] = readgeoraster(DetailedDEM);
    [DetailedREF_raster , DetailedREF_rasterref] = readgeoraster(DetailedREF);
    if idensity == 1
        [MassDensity_raster , MassDensity_rasterref] = readgeoraster(MassDensity);
    end
    
    if strcmp(DetailedDEM , TessDEM) == 0
        [TessDEM_raster , TessDEM_rasterref] = readgeoraster(TessDEM);
        [TessREF_raster , TessREF_rasterref] = readgeoraster(TessREF);
        if idensity == 1
            [TessDensity_raster, TessDensity_rasterref] = readgeoraster(TessDensity);
        end
    end
    %tiles = [2,2];
    %% write output report --- parameter definition and computation point
    fid=fopen(output_report,'w');
    p='=========================    TGF_GUI project report   =========================';
    formatSpec='%s\n';
    fprintf(fid,formatSpec,p);
    fprintf(fid,'\n');
    fprintf(fid,'%30s','Generating_data:');
    fprintf(fid,'%14s\n',datestr(now,1));
    fprintf(fid,'%30s','Product_type:');
    if ikind==1   
    fprintf(fid,'%35s\n','topographic gravitational field');
    else if ikind==2
       fprintf(fid,'%26s\n','RTM gravitational field');
    end
    end
    fprintf(fid,'%30s','Functional:');
    fprintf(fid,'%15s %d %s\n','itype=',itype,'(please refer to handbook)');
    fprintf(fid,'%30s','Gravity field type:');
    fprintf(fid,'%15s %d %s\n','ikind=',ikind,'(please refer to handbook)');
    fprintf(fid,'%30s','Earth approximation:');
    fprintf(fid,'%15s %d %s\n','iflag=',flag_earth,'(please refer to handbook)');
    fprintf(fid,'%30s','Four zones radius [deg]:');
    fprintf(fid,'%7.2f, %7.2f, %7.2f, %7.2f\n',rzones);
    fprintf(fid,'\n');
    p='----------------------  computation point information  -----------------------';
    formatSpec='%s\n';
    fprintf(fid,formatSpec,p);
    fprintf(fid,'\n');
    fprintf(fid,'%30s%10.4f\n','comlat_min:',minlat); 
    fprintf(fid,'%30s%10.4f\n','comlat_max:',maxlat); 
    fprintf(fid,'%30s%10.4f\n','comlon_min:',minlon); 
    fprintf(fid,'%30s%10.4f\n','comlon_max:',maxlon); 
    fprintf(fid,'%30s%10d\n','number of computation points:',num); 
    fprintf(fid,'\n');
    fclose(fid);      
    %%
    % Define a grid for interpolation
    figure(1)
    scatter(lon,lat,'filled')
    [lonGrid, latGrid] = meshgrid(linspace(minlon, maxlon, tiles(2)+1), linspace(maxlat, minlat, tiles(1)+1));
    figure(2)
    n=1;
    TilesIndexes = zeros(size(lat));
    t1=clock;
    for i=1:tiles(1)
        for j=1:tiles(2)
            index_points = ((lat <= latGrid(i,j)) & (lat>=latGrid(i+1,j)) &  (lon >= lonGrid(i,j)) & (lon <= lonGrid(i,j+1)));
            if sum(index_points)>=1
                fprintf('Points in tile %d : %d\n',[n,sum(index_points)])
                TilesIndexes(index_points) = n;
                
                comPointsTile = comPoints(index_points,:);
                scatter(lon(index_points),lat(index_points),'filled')
                hold on
                latlim_det = [latGrid(i+1,j) - rzones(2)*1.1 , latGrid(i,j) + rzones(2)*1.1]; 
                lonlim_det = [lonGrid(i,j) - rzones(2)*1.1 , lonGrid(i,j+1) + rzones(2)*1.1];
                %Detailed DEM/REF/density
                [DetailedDEM_crop , DetailedDEM_crop_ref] = geocrop(DetailedDEM_raster,DetailedDEM_rasterref,latlim_det,lonlim_det);
                DetailedDEM_struct = createDEMstruct(DetailedDEM_crop , DetailedDEM_crop_ref);
                [DetailedREF_crop , DetailedREF_crop_ref] = geocrop(DetailedREF_raster,DetailedREF_rasterref,latlim_det,lonlim_det);
                DetailedREF_struct = createDEMstruct(DetailedREF_crop , DetailedREF_crop_ref);
                if idensity == 1
                    [MassDensity_crop , MassDensity_crop_ref] = geocrop(MassDensity_raster , MassDensity_rasterref,latlim_det,lonlim_det);
                    MassDensity_struct = createDEMstruct(MassDensity_crop , MassDensity_crop_ref);
                    clear MassDensity_crop MassDensity_crop_ref
                else
                    MassDensity_struct = [];
                end
    
                clear DetailedDEM_crop DetailedDEM_crop_ref DetailedREF_crop DetailedREF_crop_ref
                
                % Tess DEM/REF
                latlim_tess = [latGrid(i+1,j) - rzones(3)*1.1 , latGrid(i,j) + rzones(3)*1.1]; 
                lonlim_tess = [lonGrid(i,j) - rzones(3)*1.1 , lonGrid(i,j+1) + rzones(3)*1.1];
                if strcmp(DetailedDEM , TessDEM) == 0
                    [TessDEM_crop , TessDEM_crop_ref] = geocrop(TessDEM_raster,TessDEM_rasterref,latlim_tess,lonlim_tess);
                    [TessREF_crop , TessREF_crop_ref] = geocrop(TessREF_raster,TessREF_rasterref,latlim_tess,lonlim_tess);
                    if idensity == 1
                        [TessDensity_crop , TessDensity_crop_ref] = geocrop(TessDensity_raster, TessDensity_rasterref,latlim_tess,lonlim_tess);
                    end
                else
                    [TessDEM_crop , TessDEM_crop_ref] = geocrop(DetailedDEM_raster,DetailedDEM_rasterref,latlim_tess,lonlim_tess);
                    [TessREF_crop , TessREF_crop_ref] = geocrop(DetailedREF_raster,DetailedREF_rasterref,latlim_tess,lonlim_tess);
                    if idensity ==1
                        [TessDensity_crop , TessDensity_crop_ref] = geocrop(MassDensity_raster , MassDensity_rasterref,latlim_tess,lonlim_tess);
                    end
                        
                end
                TessDEM_struct = createDEMstruct(TessDEM_crop , TessDEM_crop_ref);
                TessREF_struct = createDEMstruct(TessREF_crop , TessREF_crop_ref);
                if idensity ==1
                    TessDensity_struct = createDEMstruct(TessDensity_crop , TessDensity_crop_ref);
                    clear TessDensity_crop TessDensity_crop_ref
                else
                    TessDensity_struct = [];
                end
                clear TessDEM_crop TessDEM_crop_ref TessREF_crop TessREF_crop_ref
                [parentFolder, folderName]=fileparts(outname);
                if not (isfolder(outname))
                    mkdir(parentFolder,folderName)
                end
                outname_temp = fullfile(outname,strcat('tile_',num2str(n)));
                % Coarse DEM / REF is not clipped
                TGF_util(comPointsTile, DetailedDEM_struct, DetailedREF_struct,MassDensity_struct, TessDEM_struct, TessREF_struct, ...
                TessDensity_struct, CoarseDEM, CoarseREF, CoarseDensity,GlobalDEM,GlobalREF, outname_temp, ...
                ikind, itype,idensity, flag_earth,rzones, e);
                n = n+1;
            end
            
        end
    end
    hold off
    %length(TilesIndexes)
    %nnz(TilesIndexes)
    fprintf('Total number of non-empty tiles: %d\n',n-1)
    filename_temp = fullfile(outname,strcat('tile_',num2str(1),'.dat'));
    points_temp = load(filename_temp);
    col_size = size(points_temp,2);
    points_all = zeros(length(lat),col_size);
    index_temp = (TilesIndexes==1);
    points_all(index_temp,:)=points_temp;
    if (n-1)>1
        for i=2:n-1
            filename_temp = fullfile(outname,strcat('tile_',num2str(i),'.dat'));
            points_temp = load(filename_temp);
            index_temp = (TilesIndexes==i);
            points_all(index_temp,:)=points_temp;
            
        end
    end
    t2=clock;
    t=etime(t2,t1)/3600;
    [status, message, messageid] = rmdir(outname, 's');
    if status == 0
        disp(message)
        disp(messageid)
    end
    outname=strcat(outname,'.dat');
    fid=fopen(output_report,'a');
    p='-------------------------  integration information  --------------------------';
    fprintf(fid,'%s\n',p);
    fprintf(fid,'\n');
    fprintf(fid,'%30s%10.5f\n','Computation time[hour]:',t);
    fprintf(fid,'\n');
    p='--------------------------  Output information  -------------------------------';
    fprintf(fid,'%30s','Output file name:');
    fprintf(fid,'%15s %s \n','Output file:',outname);
    fprintf(fid,'\n');
    p='------------------------------  Results Format  ------------------------------';
    fprintf(fid,'%s\n',p);
    fprintf(fid,'\n');
    if itype == 10
        p={'lon'  'lat'  'H'  'N'  'xi'  'eta'  'dg'  'Vxx'  'Vxy'  'Vxz'  'Vyy'  'Vyz'  'Vzz'};
        formatSpec='%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s\n';
        fprintf(fid,formatSpec,p{1,:});
        p={'[deg]'  '[deg]'  '[m]'  '[m]'  '["]'  '["]'  '[mGal]'  '[E]'  '[E]'  '[E]'  '[E]'  '[E]'  '[E]'};
        formatSpec='%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s%7s\n';
        fprintf(fid,formatSpec,p{1,:});
    elseif itype == 0
        p={'lon'  'lat'  'H'  'N'};
        formatSpec='%7s%7s%7s%7s\n';
        fprintf(fid,formatSpec,p{1,:});
        p={'[deg]'  '[deg]'  '[m]'  '[m]'};
        formatSpec='%7s%7s%7s%7s\n';
        fprintf(fid,formatSpec,p{1,:});
    elseif itype==103 || itype==104 || itype==4
        p={'lon'  'lat'  'H'  'N'  'xi'  'eta'  'dg'};
        formatSpec='%7s%7s%7s%7s%7s%7s%7s\n';
        fprintf(fid,formatSpec,p{1,:});
        p={'[deg]'  '[deg]'  '[m]' '[m]'  '["]'  '["]'  '[mGal]'};
        formatSpec='%7s%7s%7s%7s%7s%7s%7s\n';
        fprintf(fid,formatSpec,p{1,:});
    elseif itype==1
        p={'lon'  'lat'  'H'  'xi'  'eta'  'dg'};
        formatSpec='%7s%7s%7s%7s%7s%7s\n';
        fprintf(fid,formatSpec,p{1,:});
        p={'[deg]'  '[deg]'  '[m]'  '["]'  '["]'  '[mGal]'};
        fprintf(fid,formatSpec,p{1,:});
    elseif itype == 2
        p={'lon'  'lat'  'H'  'xi'  'eta'  'dg'};
        formatSpec='%7s%7s%7s%7s%7s%7s\n';
        fprintf(fid,formatSpec,p{1,:});
        p={'[deg]'  '[deg]'  '[m]'  '["]'  '["]'  '[mGal]'};
        fprintf(fid,formatSpec,p{1,:});
    end
    fclose(fid);
    save(outname,'points_all','-ASCII')
end

function [struct_obj] = createDEMstruct(DEM,R)
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
    struct_obj = struct();
    struct_obj.DEM = double(flip(DEM));
    struct_obj.meta = [minlat maxlat minlon maxlon reslat reslon];
end