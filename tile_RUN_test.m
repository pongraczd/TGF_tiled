%% 

comPointsFile = "C:\PHD\2025_26_1\GeoRBF\adat\Kiskoros\topo\Mert_g_1197_topo.dat";
%comPointsFile ="C:\PHD\2025_26_1\GeoRBF\adat\grav\125_4000\faye_out2.5_JAV_H_MERIT.dat"
%DetailedDEM ="C:\PHD\2025_26_1\GeoRBF\adat\Kiskoros\topo\DEM\MERIT_DEM_Kiskoros_clipped.tif";
%DetailedREF = "C:\PHD\2025_26_1\GeoRBF\adat\Kiskoros\topo\DEM\topo_REF6500_Kiskoros_clipped.tif";
%TessDEM = "C:\PHD\2025_26_1\GeoRBF\adat\Kiskoros\topo\DEM\MERIT_DEM_Kiskoros_clipped.tif";
%TessREF = "C:\PHD\2025_26_1\GeoRBF\adat\Kiskoros\topo\DEM\topo_REF6500_Kiskoros_clipped.tif";
DetailedDEM ="C:\PHD\2025_26_1\GeoRBF\adat\MERIT_DEM_12_27_43_51.tif";
DetailedREF = "C:\PHD\2025_26_1\Grav_HUN_DATA\topo_corr\REF\TOPO_REF_6500.tif";
TessDEM = "C:\PHD\2025_26_1\GeoRBF\adat\MERIT_DEM_12_27_43_51.tif";
TessREF = "C:\PHD\2025_26_1\Grav_HUN_DATA\topo_corr\REF\TOPO_REF_6500.tif";
CoarseDEM = "C:\PHD\2025_26_1\GeoRBF\adat\MERIT_DEM_12_27_43_51_30s.tif";
CoarseREF = "C:\PHD\2025_26_1\Grav_HUN_DATA\topo_corr\REF\TOPO_REF_6500_30s.tif";
ikind = 2;
idensity = 0;
flag_earth = 0;
%% 
rzones = [0.02 0.03 0.15 0.5];
itype = 4;
e = 2.67;
outname="C:\PHD\2025_26_1\GeoRBF\adat\Kiskoros\topo\Mert_g_1197_topo_dg_corr_test";
tiles = [1,1];

TGF_tiled(comPointsFile,tiles, DetailedDEM, DetailedREF,[], TessDEM, TessREF,[], CoarseDEM, CoarseREF,[],[],[],outname,ikind, itype,idensity, flag_earth,rzones, e);
