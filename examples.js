var geometry = 
    /* color: #d63000 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[10.341854818520968, 43.6638055591542],
          [10.341854818520968, 43.514604108522185],
          [10.529995687661593, 43.514604108522185],
          [10.529995687661593, 43.6638055591542]]], null, false);

// import the library
var geeSSEBI = require('users/jpk1996jpk/geeSSEBI:geeSSEBI_library');

// First example on Landsat 8

// as input we need an input geometry (to be declared first)
// we run the appropiate function for preparing inputs 
//(there exist one specific for Landsat 5, one for Landsat 8 and one for Landsat 9)
// choose the date start and date end of the timeperiod
// the first image in this interval will be processed
var inputs_L8 = geeSSEBI.preprocessInputsL8('2020-07-01', '2020-09-01', geometry);
print(inputs_L8)
var model_result_L8 = geeSSEBI.runSSEBI(inputs_L8, geometry, 0.11, 0.11,'EPSG:32632', 100);
var daily_ET_L8 = model_result_L8.select('ET_daily');
Map.addLayer(daily_ET_L8, geeSSEBI.visParamET, 'Landsat 8 ET'); //we can add visualization parameters for customizing the view
print(daily_ET_L8);
print(model_result_L8);
var grafico_ssebi_L8 = geeSSEBI.plotSSEBI(model_result_L8, geometry, 100, 'EPSG:32632')
print(grafico_ssebi_L8);
Export.image.toDrive(model_result_L8)

var inputs_L9 = geeSSEBI.preprocessInputsL9('2023-07-01', '2023-09-01', geometry);
print(inputs_L9)
var model_result_L9 = geeSSEBI.runSSEBI(inputs_L9, geometry, 0.11, 0.11,'EPSG:32632', 100);
var daily_ET_L9 = model_result_L9.select('ET_daily');
Map.addLayer(daily_ET_L9,  geeSSEBI.visParamET, 'Landsat 9 ET'); //we can add visualization parameters for customizing the view
print(daily_ET_L9);
print(model_result_L9);
var grafico_ssebi_L9 = geeSSEBI.plotSSEBI(model_result_L9, geometry, 100, 'EPSG:32632')
print(grafico_ssebi_L9);
Export.image.toDrive(model_result_L9)
// then we run the model on the preprocessed inputs
// as arguments we specify: preprocessed inputs, geometry of interest, 
// albedo cutoff for moisture limited and radiation limited lines
// projection and spatial resolution (shall be coherent with the spatial resolution of the thermal band)



// Test the model on Landsat 5 
// The function works exactly as for Landsat 8 
var inputs_L5 = geeSSEBI.preprocessInputsL5('1991-06-01', '1991-09-01', geometry);
print(inputs_L5)
var model_result_L5 = geeSSEBI.runSSEBI(inputs_L5, geometry, 0.11, 0.11,'EPSG:32632', 100);

//We select the daily evapotranspiration band for visualizing it
var daily_ET = model_result_L5.select('ET_daily');
Map.addLayer(daily_ET, geeSSEBI.visParamET, 'Landsat 5 ET'); //we can add visualization parameters for customizing the view
print(daily_ET);
print(model_result_L5);
var grafico_ssebi_L5 = geeSSEBI.plotSSEBI(model_result_L5, geometry, 100, 'EPSG:32632')
print(grafico_ssebi_L5);
Export.image.toDrive(model_result_L5)