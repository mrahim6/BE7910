from shapely.geometry import Point
import geopandas as gpd
import pandas as pd
import rasterio
import time
# import os
# os.chdir(r"C:\Users\mrahim\Downloads\CPRA\scripts")
from helper_functions import *

#get paths
elevation_raster_path = r'C:\Users\mrahim\OneDrive - LSU AgCenter\FloodSafeHome\Risk Rating 2 Dataset\ElevationRaster'
county_shapefile_path = r'C:\Users\mrahim\OneDrive - LSU AgCenter\FloodSafeHome\Risk Rating 2 Dataset\County'
huc12_shapefile_path = r'C:\Users\mrahim\OneDrive - LSU AgCenter\FloodSafeHome\Risk Rating 2 Dataset\HUC12'
crs_shapefile_path = r'C:\Users\mrahim\OneDrive - LSU AgCenter\FloodSafeHome\Risk Rating 2 Dataset\CRS'
levee_shapefile_path = r'C:\Users\mrahim\OneDrive - LSU AgCenter\FloodSafeHome\Risk Rating 2 Dataset\Levee'
river_poly_shapefile_path = r'C:\Users\mrahim\OneDrive - LSU AgCenter\FloodSafeHome\Risk Rating 2 Dataset\RiverPoly'
flowline_shapefile_path = r'C:\Users\mrahim\OneDrive - LSU AgCenter\FloodSafeHome\Risk Rating 2 Dataset\FlowLine'
flood_depth_raster_path = r'C:\Users\mrahim\OneDrive - LSU AgCenter\FloodSafeHome\Risk Rating 2 Dataset\FloodDepth'
coastline_shapefile_path = r"C:\Users\mrahim\OneDrive - LSU AgCenter\FloodSafeHome\Risk Rating 2 Dataset\CoastLine"
building_data_path = r"C:\Users\mrahim\OneDrive - LSU AgCenter\Proposals\CPRA\Contract_RR2\Data\mp23_pdd_clara.structure_info_costs_2024.06.18 (1)"

# Load  datasets
projected_crs = 'EPSG:26915'  # Web Mercator projection (units in meters)
coastline_gdf  = gpd.read_file(coastline_shapefile_path + "/CoastLine.shp")
coastline_gdf_proj = coastline_gdf.to_crs(projected_crs)

riverpoly_gdf = gpd.read_file(river_poly_shapefile_path + "/riverPoly_FSF.shp")
riverpoly_gdf_proj = riverpoly_gdf.to_crs(projected_crs)

flowline_gdf = gpd.read_file(flowline_shapefile_path + "/FlowLine_FSF.shp")
flowline_gdf_proj = flowline_gdf.to_crs(projected_crs)

county_gdf = gpd.read_file(county_shapefile_path + "/tl_2017_LA_county.shp")
huc12_gdf = gpd.read_file(huc12_shapefile_path + "/HUC12.shp")
crs_gdf = gpd.read_file(crs_shapefile_path + "/CRS.shp")
levee_gdf = gpd.read_file(levee_shapefile_path + "/qgis/levee.shp")

elevation_dataset = rasterio.open(elevation_raster_path + "/Elevation_FSF.tif")
elevation_band = elevation_dataset.read(1)

fd_100_dataset = rasterio.open(flood_depth_raster_path + "/USFloodMap_FLUVIAL-UNDEFENDED_PERCENTILE50_2020_1in100.tif")
fd_1000_dataset = rasterio.open(flood_depth_raster_path + "/USFloodMap_FLUVIAL-UNDEFENDED_PERCENTILE50_2020_1in1000.tif")

df = pd.read_csv(building_data_path+"/mp23_pdd_clara.structure_info_costs_2024.06.18.csv")
############
start = time.time() 
geographic_attr = []
for i in range(0,len(df)):
    df1 = df.iloc[i]
    
    point = Point(df1['lon'], df1['lat'])
    point_gdf = gpd.GeoDataFrame(geometry=[point], crs='EPSG:4269')
    point_gdf_proj = point_gdf.to_crs(projected_crs)   

    county_name = get_single_value(county_gdf, point, 'NAME')
    crs = get_multiple_values(crs_gdf, point, 'CRS_Class') 
    levee = get_multiple_values(levee_gdf, point, 'systemId')
    huc12_name = get_single_value(huc12_gdf, point, 'HUC12')  
    sre = get_sre(elevation_dataset, point)
    distance_to_coast = get_distance_to_coast(point_gdf_proj, coastline_gdf_proj)   
    elevation = get_elevation(elevation_dataset, elevation_band, point)
    DTR, RiverElevation, ERR, DA, FloodDepthDifference = process_huc_data(
                                huc12_name, point, elevation, elevation_dataset, elevation_band, fd_100_dataset, fd_1000_dataset, 
                                riverpoly_gdf, flowline_gdf)

    ####
    result = {
        'structure_id': df1['structure_id'],
        'grid_point_id_2023': df1['grid_point_id_2023'],
        'geog_sre': sre,
        'geog_elev': elevation,
        'geog_countyName': county_name,
        'geog_hucName': huc12_name,
        'geog_crs': crs,
        'geog_levee_system_Id': levee,
        'geog_DTR': DTR,
        'geog_ElevationRiver': RiverElevation,
        'geog_ERR': ERR,
        'geog_DA': DA,
        'geog_FloodDepthDifference': FloodDepthDifference,
        'geog_DTC': distance_to_coast[0]
    }
    geographic_attr.append(result)

end = time.time()
print(end-start)
###############   

geographic_df = pd.DataFrame(geographic_attr)

final_df = df.merge(geographic_df, right_on = 'structure_id', right_on='structure_id')
final_df.to_csv('df_geographic_attributes.csv', index=False)


#sample script

final_df = df.merge(geographic_df, right_on = 'structure_id', right_on='structure_id')
final_df.to_csv('df_geographic_attributes.csv', index=False)

