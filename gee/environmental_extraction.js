// Load datasets
	
//WorldClim bioclimatic variables
var worldclim = ee.Image("WORLDCLIM/V1/BIO")
 .select(['bio01', 'bio02', 'bio03', 'bio04', 'bio07', 'bio12', 'bio15']);

//Tree cover (Hansen/UMD)
var treeCover = ee.Image("UMD/hansen/global_forest_change_2023_v1_11")
	.select(['treecover2000'])
	.rename('tree_cover_pct');

//Tree height (GEDI/Landsat - Potapov, 2020)
var treeHeight = ee.Image("projects/glad/GLCLU2020/Forest_height_2020")
	.rename('canopy_height_m');

//NDVI (MODIS 16-day global 250m)
var ndvi = ee.ImageCollection("MODIS/061/MOD13Q1")
 .filterDate('2023-06-01', '2023-08-31')
 .select('NDVI')
 .mean()
 .multiply(0.0001) // Scale factor for MODIS NDVI
 .rename('peak_ndvi');

//Human modification
var humanMod = ee.ImageCollection("CSP/HM/GlobalHumanModification")
 .first()
 .rename('human_modification');

// Combine datasets
var environmentalStack = worldclim
 .addBands(treeCover)
 .addBands(treeHeight)
 .addBands(ndvi)
 .addBands(humanMod);

// Extract values
var sampledObservations = environmentalStack.sampleRegions({
 collection: raccoons, 
 properties: ['id', 'latitude', 'longitude'],
 scale: 30,
geometries: true
	});
	
Export.table.toDrive({
 collection: sampledObservations,
 description: 'Raccoon_GEE_Export',
 fileFormat: 'CSV',
 selectors: [
	'id', 'bio01', 'bio02', 'bio03', 'bio04',    
	'bio07', 'bio12', 'bio15',             
	'tree_cover_pct', 'canopy_height_m', 'peak_ndvi', 
	'human_modification'      
]
});