comPointsFile = "sample_data\sample_points.dat";
DetailedDEM = "sample_data\DEM\MERIT_DEM_Kiskoros_clipped.tif";
DetailedREF = "sample_data\DEM\topo_REF6500_Kiskoros_clipped.tif";
TessDEM = "sample_data\DEM\MERIT_DEM_Kiskoros_clipped.tif";
TessREF = "sample_data\DEM\topo_REF6500_Kiskoros_clipped.tif";
CoarseDEM = "sample_data\DEM\MERIT_DEM_12_27_43_51_30s.tif";
CoarseREF = "sample_data\DEM\TOPO_REF_6500_30s.tif";
ikind = 2;
idensity = 0;
flag_earth = 0;
rzones = [0.02 0.03 0.15 0.5];
itype = 4;
e = 2.67;
outname="sample_results\sample_dg_topo";
tiles = [2,2];

TGF_tiled(comPointsFile,tiles, DetailedDEM, DetailedREF,[], TessDEM, TessREF,[], CoarseDEM, CoarseREF,[],[],[],outname,ikind, itype,idensity, flag_earth,rzones, e);
