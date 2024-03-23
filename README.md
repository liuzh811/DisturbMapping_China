# Patterns, trends, and drivers of disturbance regimes in China's forest from 1986 to 2020
## Introduction 
This repository was created to store data and codes used for Liu et al., 2023, Communications Earth & Environment.
## Citation
Liu, Z., Wang, W.J., Ballantyne, A. et al. Forest disturbance decreased in China from 1986 to 2020 despite regional variations. Commun Earth Environ 4, 15 (2023). https://doi.org/10.1038/s43247-023-00676-x
## folder description
1. Data: contains data to generate the figures in the manuscript.
2. Code: contains R codes for analysis and plot.
3. Disturbance_layer: contains annually % disturbed forest at 250m resolution. If you need 30-m resolution data, contact Z. Liu.
## Abstract
Human activities have altered disturbance patterns in many parts of world, but there is no quantitative information on patterns and trends of forest disturbance regimes in China. Here we applied a spectral-temporal segmentation approach over all available Landsat data to map individual disturbance patches and characterize the patterns and trends in disturbance rate, size, frequency, and severity across China’s forests. From 1986 to 2020, about 39.7% of China’s forests were disturbed with an annual rate of 1.16 ± 0.41% yr-1. The disturbance decreased at a rate of -390 ± 142 km2 yr-1, primarily driven by the effective implementation of forest protection policy since 2000s. The rate, frequency, and size of disturbance generally intensified in Southeast, but weakened in Northeast China. Our high-quality, spatially explicit disturbance map provides an essential data layer to understand the landscape-scale drivers of forest dynamics and functions for important but less understood pan-temperate forest regions. 

## Access to 30m disturbance data
### Load data
var c01 = ee.Image("users/liuzh822/DisturbMap/FinalClass_r01");

var c02 = ee.Image("users/liuzh822/DisturbMap/FinalClass_r02");

var c03 = ee.Image("users/liuzh822/DisturbMap/FinalClass_r03");

var c04 = ee.Image("users/liuzh822/DisturbMap/FinalClass_r04");

var c05 = ee.Image("users/liuzh822/DisturbMap/FinalClass_r05");

var c06 = ee.Image("users/liuzh822/DisturbMap/FinalClass_r07");

var c07 = ee.Image("users/liuzh822/DisturbMap/FinalClass_r08");

var c08 = ee.Image("users/liuzh822/DisturbMap/FinalClass_r7");

var c09 = ee.Image("users/liuzh833/DisturbMap/FinalClass_r10");

### Mosaic
var classes = ee.ImageCollection([c01, c02, c03, c04, c05, c06,c07,c08,c09]);

var classes = classes.mosaic().selfMask().rename("lossyear");
 
/*

value 1-36: fire disturbance 1985-2020; 

101-136: forest management-related disturbance 1985-2020; 

201-236: other disturbance 1985-2020

*/


var fire_dis = classes.updateMask(classes.lte(36).and(classes.gte(1))).add(1984);

var forest_dis = classes.updateMask(classes.lte(136).and(classes.gte(101))).add(-100).add(1984);

var other_dis = classes.updateMask(classes.lte(230).and(classes.gte(201))).add(-200).add(1984);

var all_dis = fire_dis.unmask().add(forest_dis.unmask()).add(other_dis.unmask()).clip(chn).selfMask();


### Visualization

var startYear = 1985;

var endYear = 2020;

// red --> yellow --> blue [1985 --> 2020]

var palette = ['#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'];

var yodVizParms = {

  min: startYear,
  
  max: endYear,
  
  palette: palette
  
};


Map.addLayer(fire_dis.selfMask(), yodVizParms, 'fire disturbance', false);

Map.addLayer(forest_dis.selfMask(), yodVizParms, 'forest management disturbance', false);

Map.addLayer(other_dis.selfMask(), yodVizParms, 'other disturbance', false);

Map.addLayer(all_dis.selfMask(), yodVizParms, 'all disturbance', false);


