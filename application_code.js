var geometry = /* color: #000000 */ee.Geometry.MultiPoint();
var bckcolor = "#FFFFFF"
var ERA5H = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY");
var ERA5D = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR");
var scale = 100;
//var visParamET = {bands : ['ET_daily'], max : 5.5, min: 0, opacity :1, palette: ['0af4ff', '39ff51', '98ff4b', 'ffe94d', 'ff7f23', 'ff1004']};
//var visParamNDVI = {bands : ['NDVI'], max : 1, min: 0, opacity :1, palette: ['#DEB887', '#BDB76B', '#8FBC8F', '#556B2F', '#006400']};
var visParamET = {"opacity":1,"bands":["ET_daily"],"min":0,"max":7,"palette":['0af4ff', '39ff51', '98ff4b', 'ffe94d', 'ff7f23', 'ff1004']};
var visParamNDVI = {bands : ['NDVI'], max : 1, min: 0, opacity :1, palette: ['#DEB887', '#BDB76B', '#8FBC8F', '#556B2F', '#006400']};
var visParamGeo = {fillColor : '00000000'};

// First visualization settings
Map.setOptions('HYBRID');
Map.setCenter(12.4964, 41.9028, 8); // Nos centramos en Roma, Italia
Map.setControlVisibility(false, true, true, true, true, false, false);


/*  APP STRUCTURE
-addition of drawing tools
-declaration of functions for running SSEBI
-SSEBI is run immediately on draw
*/

/* ========= Function definition section ============
Hereafter, the base functions are defined. 
1 - maskL8sr Masks the image, applies the scaling factors for reflectance.
2 - compute NDWI computes the NDWI, for later applying the masking (removal)
    of water bodies (NDWI > 0.2)
3 - maskWater removes the water bodies. 
4 - getDateProperty gets the date, for taking the ERA5 data of the same day 

Later, the higher level functions are defined. 
1 - preprocessInputs
2 - runSSEBI

preprocessInputs performs the following steps:
1) declaration of functions for NDVI, fraction cover, emissivity, albedo computation.
2) LST is used to compute the instantaneous emitted radiation.
3) incoming radiation (both instantaneous and daily) is retrieved from ERA5
4) The net Radiation is computed for both timesteps. 
5) Only the useful bands are maintained

runSSEBI runs the rest of the model, computing the slope of 
the radiation limited, and moisture limited lines. Then computes Twet and Thot for each pixel, 
from its albedo, and finally the evaporative fraction.
Moreover, soil Flux coefficient is computed from fraction cover
And later used to derive the soil flux from radiation. 
The remaining energy of radiation is partitioned according to the evaporative fraction.
Instantaneous ET is computed by dividing by the latent heat flux.
Hourly ET is obtained by multiplying the instantaneous one by 3600.
Daily ET is obtained by multiplying the instantaneous one by the ratio of daily to instantaneous 
Radiation (CDI).
Latent heat of evapotranspiration is assumed constant (2.46 MJ/Kg)

For simplicity, the final output contains the bands:
1 - albedo
2 - surface temperature (K)
3 - emissivity
4 - NDVI
5 - fraction cover
6 - Net Radiation Instantaneous (W/m^2)
7 - LE instantaneous (W/m^2)
8 - H instantanous (W/m^2)
9 - G (soil flux) instantaneous (W/m^2)
10 - Daily ET (mm/day)
11 - Timestamp
*/



var maskL8sr = function(image) {
    // Bit 0 - Fill
    // Bit 1 - Dilated Cloud
    // Bit 2 - Cirrus
    // Bit 3 - Cloud
    // Bit 4 - Cloud Shadow
    var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
    var saturationMask = image.select('QA_RADSAT').eq(0);
    
    //My own cloud distance mask to be more restrictive
    var cdist = image.select("ST_CDIST").multiply(0.01);
    var cdist_mask = cdist.gt(ee.Image.constant(10))
    
    // Apply the scaling factors to the appropriate bands.
    var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
    var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
    var constant = ee.Image.constant(0.0001).rename('area')
    // Aggiungere maschera per togliere valori di ST eccessivi, tipo maggiori di 350 K e minori di 273?
    // Replace the original bands with the scaled ones and apply the masks.
    var im_masked = image.addBands(opticalBands, null, true)
      .addBands(thermalBands, null, true)
      .addBands(constant)
      .updateMask(qaMask)
      .updateMask(saturationMask);
    return im_masked.updateMask(im_masked.select('ST_B10').gt(ee.Image.constant(273)))
      .updateMask(im_masked.select('ST_B10').lt(ee.Image.constant(340)));
  };
  
var computeNDWI = function(image){
  var NDWI = image.select("SR_B3").subtract(image.select("SR_B6"))
  .divide(image.select("SR_B3").add(image.select("SR_B6"))).rename("NDWI")
  return image.addBands(NDWI)
}

var maskWater = function(image){
  return image.updateMask(image.select("NDWI").lt(0.0))
}
    
var getDateProperty = function(feature) {
  return feature.get('date');
};

// Preprocess inputs function

var prepareInputs = function (startDate,
    endDate,
    geometry,
    collection) {
    // declare constant
    var B_CONST = ee.Number(5.67).multiply(ee.Number(10).pow(-8));
    //import L8
    var L8 = ee.ImageCollection(collection).filterBounds(geometry)
    // Delcare the cloud masking function and applying scaling factors to the reflectance values
    
    // removed maskL8sr duplicate
    // put again if needed

    // Declare functions for computing NDVI and other spectral properties
    // They are normalized to be restricted into the range between 0 and 1

    function obtainNDVI(image) {
        var NDVI = image.select('SR_B5').subtract(image.select('SR_B4'))
            .divide(image.select('SR_B5').add(image.select('SR_B4')))
            .rename('NDVI');
        var NDVI_ok = NDVI.multiply(NDVI.gt(ee.Image.constant(0))).multiply(NDVI.lt(ee.Image.constant(1)))
            .add(NDVI.gt(ee.Image.constant(1))).rename("NDVI")
        var result = image.addBands(NDVI_ok);
        return (result)
    }


    function obtainFC(image) {
        var FC = image.select("NDVI").subtract(ee.Image.constant(0.2))
            .divide(ee.Image.constant(0.6)).pow(ee.Image.constant(2)).rename("FC");
        FC = FC.where(FC.gt(1), ee.Image.constant(1)).where(image.select("NDVI").lt(0.2), ee.Image.constant(0))
        var result = image.addBands(FC);
        return (result)
    }

    function obtainALBEDO(image) {
        var ALBEDO = image.select('SR_B1').multiply(ee.Image.constant(0.13))
            .add(image.select('SR_B2').multiply(ee.Image.constant(0.115)))
            .add(image.select('SR_B3').multiply(ee.Image.constant(0.143)))
            .add(image.select('SR_B4').multiply(ee.Image.constant(0.18)))
            .add(image.select('SR_B5').multiply(ee.Image.constant(0.281)))
            .rename('albedo');
        var ALBEDO_ok = ALBEDO.multiply(ALBEDO.gt(ee.Image.constant(0))).multiply(ALBEDO.lt(ee.Image.constant(1)))
            .add(ALBEDO.gt(ee.Image.constant(1)))
        var result = image.addBands(ALBEDO_ok.rename('albedo'));
        return (result)
    }
    //rewrite and add
    function obtainEMISS(image) {
        var emiss = ee.Image.constant(0.971).multiply(ee.Image.constant(1).subtract(image.select("FC")))
            .add(ee.Image.constant(0.982).multiply(image.select("FC"))).rename("emiss");
        var result = image.addBands(emiss);
        return (result);
    }
    /* ========== Preprocess function ============
    Produces the albedo discrete band, for the computation of the dry and wet regression lines
    Renames bands and selects the bands necessary for running the SEBS model.
    */
    var preprocess = function (image) {
        var albedo_disc = image.select('albedo').multiply(ee.Image.constant(100)).round().rename('albedo_disc');
        var result = image.select('albedo')
            .addBands(image.select("NDVI"))
            .addBands(image.select("FC"))
            .addBands(image.select('emiss'))
            .addBands(image.select('ST_B10').rename('LST'))
            .addBands(albedo_disc);
        return (result);
    };

    // function for emitted radiation
    var computeEmitted = function (image) {
        var emitted = image.select('LST').pow(4);
        emitted = emitted.multiply(image.select('emiss')).multiply(ee.Image.constant(B_CONST)).rename("longwaveEmitted");
        var result = image.addBands(emitted);
        return (result);
    };


    /* 
    a function to obtain the radiation will be necessary
    ERA-5 might be a suitable universal source
    */

    var shortwaveNetI = function (image) {
        var result = image.select("ShortwaveDownI").multiply(ee.Image.constant(1).subtract(image.select("albedo"))).rename("ShortwaveNetI")
        return (image.addBands(result))
    }
    var longwaveNetI = function (image) {
        var pt1 = image.select("LongwaveDownI").multiply(image.select('emiss')).subtract(image.select('longwaveEmitted')).rename("LongwaveNetI")//.multiply(image.select("emiss")).rename("LI")
        //var pt2 = image.select("emiss").multiply(image.select("LST").pow(ee.Image.constant(4))).multiply(ee.Image.constant(B_CONST)).rename("LO")
        //var mean_emitted = pt2.reduceRegion(ee.Reducer.mean(),geometry, 100, "EPSG:32632")
        //print("Longwave emitted")
        //print(mean_emitted.get("LO"))
        //var result = pt1.select("LI").subtract(pt2.select("LO")).rename("LongwaveNetI")
        return (image.addBands(pt1))
    }
    var shortwaveNetD = function (image) {
        var result = image.select("ShortwaveDownD").multiply(ee.Image.constant(1).subtract(image.select("albedo"))).rename("ShortwaveNetD")
        return (image.addBands(result))
    }
    var longwaveNetD = function (image) {
        var pt1 = image.select("LongwaveDownD").rename("LongwaveNetD")//.multiply(image.select("emiss"))
        //.rename("LI")
        //var pt2 = image.select("emiss").multiply(image.select("LST").pow(ee.Image.constant(4))).multiply(ee.Image.constant(B_CONST)).rename("LO")
        //var result = pt1.select("LI").subtract(pt2.select("LO")).rename("LongwaveNetD")
        return (image.addBands(pt1))
    }

    var computeNetI = function (image) {
        var result = shortwaveNetI(image).select("ShortwaveNetI").add(longwaveNetI(image).select("LongwaveNetI"))
            .rename("RnI")
        return (image.addBands(result))
    }

    var computeNetD = function (image) {
        var result = shortwaveNetD(image).select("ShortwaveNetD").multiply(ee.Image.constant(1).subtract(image.select('albedo'))).add(longwaveNetD(image).select("LongwaveNetD"))
            .rename("RnD")
        return (image.addBands(result))
    }

    var obtainRadiationI = function (image,
        startDate,
        endDate,
        hour) {
        //get the radiation of the day and hour corresponding to image acquisition
        var targetRN = ERA5H.filterDate(image.date(), image.date().advance(1, "day")).filter(ee.Filter.eq("hour", hour)).first();

        var short = targetRN.select("surface_solar_radiation_downwards_hourly")
            .divide(ee.Image.constant(3600)).rename("ShortwaveDownI");// dividi per 3600 per ottenere istantanea

        var longwave = targetRN.select("surface_thermal_radiation_downwards_hourly").divide(ee.Image.constant(3600)).rename("LongwaveDownI");
        var result = computeNetI(image.addBands(short.clip(geometry)).addBands(longwave.clip(geometry)));
        return (result)
    }
    var obtainRadiationD = function (image,
        startDate,
        endDate) {
        // get the radiation corresponding to the day of the image
        var targetRN = ERA5D.filterDate(image.date(), image.date().advance(1, "day")).first();
        print("daily ERA5")
        var short = targetRN.select("surface_solar_radiation_downwards_sum")
            .rename("ShortwaveDownD");
        var longwave = targetRN.select("surface_net_thermal_radiation_sum").rename("LongwaveDownD");
        var result = computeNetD(image.addBands(short.clip(geometry)).addBands(longwave.clip(geometry)));
        return (result)
    }

    // final step where we apply the functions declared above
    var an_image = L8.filterDate(startDate, endDate)
        .filterBounds(geometry).map(maskL8sr)
        .map(computeNDWI).map(maskWater)
        .map(obtainALBEDO)
        .map(obtainNDVI).map(obtainFC).map(obtainEMISS)
        .map(preprocess).map(computeEmitted)
        .first().clip(geometry);

    // hour computation
    var time = an_image.date()// get the time
    var hour = time.millis().divide(60 * 60 * 1000 * 24) //convert time to millis
        .subtract(time.millis().divide(60 * 60 * 1000 * 24).round()).multiply(24).round() //subtract the day to get the hour
    // finally apply the functions declared above
    var result = obtainRadiationD(obtainRadiationI(an_image,
        startDate,
        endDate,
        hour),
        startDate,
        endDate)
  var title = ui.Label('Image acquired on ' + an_image.get('DATE_ACQUIRED').getInfo(), {fontSize: '20px', color: 'black', backgroundColor: 'rgba(255, 255, 255, 0.8)', padding: '10px 16px'});
  Map.add(title);


    return (result)

}

// Main function Pt. 2: declaration of runSSEBI
function runSSEBI(image,
    geometry,
    rt,
    mt,
    CRS,
    scale) {
    // get the projection of the image for the subsequent computations 
    image = image.reproject({crs: CRS,
        scale : scale
      })                 
    /* Soil flux formulation according to Su 2002
    G0 = Rn · [Γc + (1 − fc) · (Γs − Γc)]
    (Γc = 0.05)
    (Γs = 0.315)
    */
    var addSFCoef = function (image) {
        var sfCoef = ee.Image.constant(0.05).add(ee.Image.constant(0.265).multiply(ee.Image.constant(1).subtract(image.select('FC'))));
        var result = image.addBands(sfCoef.rename('sfCoef'));
        return (result);
    }
    // declare the function for the final computation of ET
    var computeET = function (image) {
        var soilFlux = image.select('RnI').multiply(image.select('sfCoef'));
        var remainingRn = image.select('RnI').subtract(soilFlux)
        var H = remainingRn.multiply(ee.Image.constant(1).subtract(image.select('EF')));
        var LE = remainingRn.multiply(image.select('EF'));
        var CDI = image.select('ShortwaveDownD').divide(image.select("ShortwaveDownI")).rename("CDI"); // conversion factor of instantaneous to daily
        //var CDI2 = image.select('ShortwaveDownD').divide(image.select("ShortwaveDownI")).rename("CDI2")
        // Rivedere calore latente di evaporazione ma dovrebbe essere Ok
        var ET_Kg_sm2 = LE.divide(ee.Image.constant(2464705));//value in J to yield millimetersmeters per second unit mm/s
        var ET_mm_h = ET_Kg_sm2.multiply(ee.Image.constant(3600));//convert to mm/h
        var ET_daily = ET_Kg_sm2.multiply(CDI).rename("ET_daily"); //Convert instant to daily following Sobrino
        var result = image.addBands(soilFlux.rename('soilFlux')) // add all the bands of the newly computed values
            .addBands(H.rename('H')).addBands(LE.rename('LE')).addBands(ET_Kg_sm2.rename('ET_Kg_sm2'))
            .addBands(ET_mm_h.rename('ET_mm_h')).addBands(CDI).addBands(ET_daily);
        return (result);
    }
    // determine the evaporative fraction, by linear regression of the minimum and maximum LST

    var computeEF = function (image) {

        var hist = image.select('albedo_disc').reduceRegion({
            reducer: ee.Reducer.frequencyHistogram(),
            geometry: geometry,
            scale: scale,
            crs: CRS
        }
        )
        print(hist.get('albedo_disc'))
        var classes = ee.Dictionary(hist.get('albedo_disc')).keys()

        var maskDiscrete = function (value, image) {
            image = ee.ImageCollection(image)
            var mask = ee.Image(image.first().select('albedo_disc')).eq(value).rename(ee.String('EQ').cat(ee.Number(value).format('%.0f')));
            var result = image.first().mask(mask).reduceRegion(
                {
                    reducer: ee.Reducer.min(),
                    geometry: geometry,
                    scale: scale

                })
            return (result)
        }


        var create_val_bands = function (image) {
            var hist = image.select('albedo_disc').reduceRegion({
                reducer: ee.Reducer.frequencyHistogram(),
                geometry: geometry,
                scale: scale,
                crs: CRS
            }
            )
            var conv2n = function (string) {
                var result = ee.Number.parse(string)
                return (result)
            }
            var classes = ee.Dictionary(hist.get('albedo_disc')).keys()
            var numeric_classes = classes.map(conv2n)
            var final = numeric_classes.map(function (value) {
                var result = image.addBands(image.select('albedo_disc').eq(ee.Image.constant(ee.Number.parse(value).int())).rename(ee.String('EQ').cat(ee.Number(value).format('%.0f')))).set('class', ee.Number(value).divide(100));
                return (result)
            })


            var final_result = ee.ImageCollection(final).map(function (image) {
                var maschera = image.select('EQ.*');
                var masked = image.mask(maschera);

                var albedoMin = masked.select('albedo').reduceRegion({
                    reducer: ee.Reducer.min(),// can be replace
                    geometry: geometry,
                    scale: scale,
                    maxPixels: 1e10
                })
                var albedoMax = masked.select('albedo').reduceRegion({
                    reducer: ee.Reducer.max(),
                    geometry: geometry,
                    scale: scale,
                    maxPixels: 1e10
                })
                var lstMin = masked.select('LST').reduceRegion({
                    reducer: ee.Reducer.min(),
                    geometry: geometry,
                    scale: 100,
                    maxPixels: 1e10
                })
                var lstMax = masked.select('LST').reduceRegion({
                    reducer: ee.Reducer.max(),
                    geometry: geometry,
                    scale: scale,
                    maxPixels: 1e10
                })


                albedoMin = ee.List([albedoMin.get('albedo'), -9999])
                    .reduce(ee.Reducer.firstNonNull());
                albedoMax = ee.List([albedoMax.get('albedo'), -9999])
                    .reduce(ee.Reducer.firstNonNull());
                lstMin = ee.List([lstMin.get('LST'), -9999])
                    .reduce(ee.Reducer.firstNonNull());
                lstMax = ee.List([lstMax.get('LST'), -9999])
                    .reduce(ee.Reducer.firstNonNull());


                // Create a feature with null geometry and NDVI value and date as properties
                var f = ee.Feature(null, {
                    'class': image.get('class'),
                    'albedoMin': albedoMin,
                    'albedoMax': albedoMax,
                    'lstMin': lstMin,
                    'lstMax': lstMax,
                    'constant': 1
                })
                return f
            })
            return (final_result)
        }

        var prova = create_val_bands(image)

        var lowLR = ee.Dictionary(prova.filter(ee.Filter.gt('albedoMin', rt)).reduceColumns({
            reducer: ee.Reducer.linearRegression({
                numX: 2,
                numY: 1
            }),
            selectors: ['constant', 'albedoMin', 'lstMin']
        }));

        var highLR = ee.Dictionary(prova.filter(ee.Filter.gt('albedoMax', mt)).reduceColumns({
            reducer: ee.Reducer.linearRegression({
                numX: 2,
                numY: 1
            }),
            selectors: ['constant', 'albedoMax', 'lstMax']
        }));
        // High is the upper line (max sensible Heat)
        //Low is the lower line (max latent Heat)
        var coefLR = ee.Array(lowLR.get('coefficients')).toList();
        var coefHR = ee.Array(highLR.get('coefficients')).toList();
        // Extract the y-intercept and slope.
        var yIntLR = ee.List(coefLR.get(0)).get(0); // y-intercept
        var slopeLR = ee.List(coefLR.get(1)).get(0); // slope

        var yIntHR = ee.List(coefHR.get(0)).get(0); // y-intercept
        var slopeHR = ee.List(coefHR.get(1)).get(0); // slope

        var EF = image.select('albedo').multiply(ee.Image.constant(slopeHR)).add(ee.Image.constant(yIntHR)).subtract(image.select('LST'))
            .divide(image.select('albedo').multiply(ee.Image.constant(slopeHR)).add(ee.Image.constant(yIntHR)).subtract(image.select('albedo').multiply(ee.Image.constant(slopeLR)).add(ee.Image.constant(yIntLR)))).rename('EF')
        var proiezione = EF.select('EF').projection().getInfo();
        var result = image.addBands(EF.rename('EF'))
            .addBands(ee.Image.constant(ee.Number(yIntLR)).rename('yIntLR').reproject(proiezione.crs))
            .addBands(ee.Image.constant(ee.Number(slopeLR)).rename('slopeLR').reproject(proiezione.crs))
            .addBands(ee.Image.constant(ee.Number(yIntHR)).rename('yIntHR').reproject(proiezione.crs))
            .addBands(ee.Image.constant(ee.Number(slopeHR)).rename('slopeHR').reproject(proiezione.crs));
        return (result)
    }
    // Finally compute
    var provaET = computeET(addSFCoef(computeEF(image)));
    var proj = provaET.select('NDVI').projection().getInfo();
    provaET = provaET.addBands(ee.Image.constant(provaET.date().millis().divide(1000)).rename('date_seconds').reproject(proj.crs))
    
    return (provaET)
}
// Function for final output selection
var selectOUT = function(image){
  return image.select('albedo', 'LST','emiss', 'NDVI', 'FC', 'RnI', 'LE', 'H', 'soilFlux', 'ET_daily', 'date_seconds')
}

// Plotting functions, for generating the nice SSEBI plot, and histogram of ET.
var plotSSEBI = function(image, rg, scale, prj){
    
var SAMPLES = image.sample(
    {region : rg, 
    scale : scale , 
    projection: prj, 
    numPixels: 1000, 
    seed : 2}).map(function(f){
      f = f.set('group', 'real data')
      return f
    })
  //print(SAMPLES)
  
  
  var MOISTURE = image.addBands(image.select('albedo')
            .multiply(image.select('slopeHR'))
            .add(image.select('yIntHR')).rename('LST'),  null,  true)
  .sample(
    {region : rg, 
    scale : scale , 
    projection: prj, 
    numPixels : 1000, 
    seed : 2}).map(function(f){
      f = f.set('group', 'moisture')
      return f
    })
  
  
  var RADIATION = image.addBands(image.select('albedo').multiply(image.select('slopeLR')).add(image.select('yIntLR')).rename('LST'),
  null,  true)
  .sample(
    {region : rg, 
    scale : scale , 
    projection: prj, 
    numPixels : 1000, 
    seed : 2}).map(function(f){
      f = f.set('group', 'radiation')
      return f
    })
  
  var data4chart = SAMPLES.merge(RADIATION).merge(MOISTURE)
  
  
  
  var ssebichart = ui.Chart.feature.groups(data4chart, 'albedo', 'LST', 'group')
  .setChartType('ScatterChart')
    .setOptions({title: 'SSEBI plot',
                    hAxis: {
                      title: 'Albedo',
                      titleTextStyle: {italic: false, bold: true}
                    },
                    vAxis: {
                      title: 'LST (K)',
                      titleTextStyle: {italic: false, bold: true}
                    },
                    series: {
                      0 : {pointSize: 0.1},
                    1 : {pointSize : 0,
                    lineWidth : 2},
                    2 : {pointSize : 0,
                    lineWidth : 2}}
                  })
  return(ssebichart)
}

var etHist = function(image, rg, scale){
    var et_chart = ui.Chart.image.histogram({
        image: image.select('ET_daily'), 
        region: rg, 
        scale: scale, 
        maxPixels: 10000000000000}
        ).setOptions({title: 'ET histogram',
            hAxis: {
              title: 'ET (mm/day)',
              titleTextStyle: {italic: false, bold: true}
            },
            vAxis: {
              title: 'Frequency',
              titleTextStyle: {italic: false, bold: true}
            }
          })
    return(et_chart)
}

// var clarificatio
// Label for "Bands in output:" with bold and centered style
var titleLabel = ui.Label({
  value: 'Bands in output:',
  style: {
    fontSize: '14px', 
    fontWeight: 'bold', 
    textAlign: 'center',
    margin: '1px 0 1px 140px'  // Adding some space at the bottom
  }
});

// Label for the rest of the content
var contentLabel = ui.Label({
  value: '1 - Albedo\n2 - Surface Temperature (K)\n3 - Emissivity\n4 - NDVI\n5 - Fraction Cover\n6 - Net Radiation Instantaneous (W/m^2)\n7 - LE instantaneous (W/m^2)\n8 - H instantaneous (W/m^2)\n9 - G (soil flux) instantaneous (W/m^2)\n10 - Daily ET (mm/day)\n11 - Acquisition time (seconds)',
  style: {
    fontSize: '12px', 
    height: '160px', 
    whiteSpace: 'pre-wrap'
  }
});

// Panel to contain both labels
var clarificationPanel = ui.Panel({
  widgets: [titleLabel, contentLabel], 
  style: {width: '250px'}  // Adjust width as necessary
});


var loadInputs = function () {

    var startDate = startDateTexbox.getValue();
    var endDate = new Date();
    var collection = collectionSelector.getValue();
    // geometry is taken in a different way
    var rt = RadLimSlider.getValue();
    var mt = MoistLimSlider.getValue();


    return {
        startDate: startDate,
        endDate: endDate,
        collection: collection,
        rt: rt,
        mt: mt
    };
};


/*
// Pannello di errore creato separatamente
var errorPanel = ui.Panel({
    style: {
        position: 'bottom-center',
        padding: '8px',
        backgroundColor: '#FFFFFF',
        border: '1px solid red'
    },
    widgets: [
        ui.Label({
            value: "ERROR: NO GEOMETRY DRAWN. PLEASE DRAW A NEW GEOMETRY.",
            style: {
                fontSize: '15px',
                fontWeight: 'bold',
                color: 'red',
                padding: '10px',
                margin: '0px'
            }
        })
    ]
});*/

// Pannello di errore creato separatamente con messaggio più scorrevole
var errorPanel = ui.Panel({
    style: {
        position: 'top-center', // Cambia la posizione in alto
        padding: '10px',
        backgroundColor: '#f8d7da', // Colore di sfondo più tenue
        border: '0.5px solid #e07a7f', // Bordo in linea con il colore dello sfondo
        width: '300px'
    },
    widgets: [
        ui.Label({
            value: "Oops! No geometry detected. Please draw a new area on the map to proceed.",
            style: {
                fontSize: '16px', // Font size leggermente ridotto per scorrevolezza
                fontWeight: 'normal', // Peso font normale per rendere la lettura più fluida
                color: '#721c24', // Colore rosso scuro per visibilità
                padding: '10px',
                textAlign: 'center', // Allinea il testo al centro
                margin: '0px'
            }
        })
    ]
});


// Funzione per controllare se il poligono è vuoto
var checkEmptyPolygon = function(polygon) {
    var polygonCoordinates = polygon.coordinates().getInfo(); // Ottiene le coordinate
    var numVertices = polygonCoordinates.length; // Conta i vertici
    
    if (numVertices === 0) {
        return errorPanel; // Restituisce il pannello di errore se il poligono è vuoto
    } else {
        return null; // Se ci sono vertici, non serve l'errore
    }
};


var runModel = function () {
  
    Map.clear();
    Map.setOptions('HYBRID');
  
    // Verifica se esiste un'area di interesse
    var aoiLayer = drawingTools.layers().get(0); 
    // Ottieni il poligono (AOI) e verifica se è vuoto
    var polygon = aoiLayer.getEeObject();
    if (checkEmptyPolygon(polygon) !== null) {
        // Se il poligono è vuoto, mostra l'errore e termina
        Map.add(checkEmptyPolygon(polygon));
        return;
    }


    // Make the chart panel visible the first time a geometry is drawn.
    if (!chartPanel.style().get('shown')) {
        chartPanel.style().set('shown', true);
    }
  
    var waitingMessage = ui.Label('Model is running...')
    chartPanel.widgets().set(1, waitingMessage)
    
    // Get the drawn geometry; it will define the reduction region.
    var aoi = drawingTools.layers().get(0).getEeObject();
    Map.centerObject(aoi);
    var inputs = loadInputs()
    // Set the drawing mode back to null; turns drawing off.
    drawingTools.setShape(null);
    var inp = prepareInputs(inputs.startDate, inputs.endDate, aoi, inputs.collection); // fix to 
    var projection = inp.select('NDVI').projection().getInfo();
    print(projection)
    var ssebi_out = runSSEBI(inp, aoi, inputs.rt, inputs.mt, projection.crs, 100);
    print(ssebi_out);
    // add the plots to the plot panel
    
    
    var SSEBI_plot = plotSSEBI(ssebi_out, aoi, 100, projection.crs)
    var hist_ET = etHist(ssebi_out, aoi, 100)
    
    // Download button
    var filename = 'Result_image';
    var url = ee.Image(selectOUT(ssebi_out)).getDownloadURL({
      'region': aoi,
      'scale': 100,
      'crs': projection.crs,
      'name': filename,
      'filePerBand': false, 
      'format': 'ZIPPED_GEO_TIFF'
    });
    var downloadPanel = ui.Panel({style: {stretch: 'horizontal', position: 'bottom-center'}});
    downloadPanel.add(ui.Label(
        'DOWNLOAD THE PRODUCT', 
        {},url
    ));
    Map.add(downloadPanel);
    



// Function to create a legend for a specific index
var paleta_ET = ['0af4ff', '39ff51', '98ff4b', 'ffe94d', 'ff7f23', 'ff1004'];
// Creates a color bar thumbnail image for use in legend from the given color palette
function makeColorBarParams(palet)
 {
  return {
    bbox: [0, 0, 1, 0.1],
    dimensions: '100x10',
    format: 'png',
    min: 0,
    max: 1,
    palette: palet//['0af4ff','39ff51','98ff4b','ffe94d','ff1004','ff7f23'],
  };
}

  
// Create a panel with three numbers for the legend
    var legendLabels_et = ui.Panel({
      widgets: [
        ui.Label(visParamET.min, {margin: '4px 8px'}),
        ui.Label(
            ((visParamET.max-visParamET.min) / 2+visParamET.min),
            {margin: '4px 8px', textAlign: 'center', stretch: 'horizontal'}),
        ui.Label(visParamET.max, {margin: '4px 8px'})
      ],
      layout: ui.Panel.Layout.flow('horizontal')
    });
    
    // Legend title
    var legendTitle_et = ui.Label({
      value: 'ET daily (mm/day)',
      style: {fontWeight: 'bold',  fontSize: '15px',
    margin: 'auto',
    padding: '4px'}
    });
    
    // Create the colour bar for the ndvi legend
    var colorBar_et = ui.Thumbnail({
      image: ee.Image.pixelLonLat().select(0),
      params: makeColorBarParams(paleta_ET),
      style: {stretch: 'horizontal', margin: '0px 8px', maxHeight: '24px'},
    });
    
    // Add the legendPanel to the map
    var legendPanel_et = ui.Panel([legendTitle_et, colorBar_et, legendLabels_et]);
    legendPanel_et.style().set({
     position: 'bottom-right', // Posiziona in basso a destra
     padding: '8px',
     width: '340px'
    });







var paleta_NDVI = ['#DEB887', '#BDB76B', '#8FBC8F', '#556B2F', '#006400'];
var visParamNDVI = {bands : ['NDVI'], max : 1, min: 0, range: [-1, 1], opacity :1, palette: ['#DEB887', '#BDB76B', '#8FBC8F', '#556B2F', '#006400']};

  
// Create a panel with three numbers for the legend
    var legendLabels_ndvi = ui.Panel({
      widgets: [
        ui.Label(visParamNDVI.min, {margin: '4px 8px'}),
        ui.Label(
            ((visParamNDVI.max-visParamNDVI.min) / 2+visParamNDVI.min),
            {margin: '4px 8px', textAlign: 'center', stretch: 'horizontal'}),
        ui.Label(visParamNDVI.max, {margin: '4px 8px'})
      ],
      layout: ui.Panel.Layout.flow('horizontal')
    });
    
    // Legend title
    var legendTitle_ndvi = ui.Label({
      value: 'NDVI',
      style: {fontWeight: 'bold',  fontSize: '15px',
    margin: 'auto',
    padding: '4px'}
    });
    
    // Create the colour bar for the ndvi legend
    var colorBar_ndvi = ui.Thumbnail({
      image: ee.Image.pixelLonLat().select(0),
      params: makeColorBarParams(paleta_NDVI),
      style: {stretch: 'horizontal', margin: '0px 8px', maxHeight: '24px'},
    });
    
    // Add the legendPanel to the map
    var legendPanel_ndvi = ui.Panel([legendTitle_ndvi, colorBar_ndvi, legendLabels_ndvi]);
    legendPanel_ndvi.style().set({
     position: 'bottom-right', // Posiziona in basso a destra
     padding: '8px',
     width: '340px'
    });





    

    // Add the plots and the button to the output panel
    chartPanel.clear()
    chartPanel.widgets().set(2, clarificationPanel);
    chartPanel.widgets().set(3, SSEBI_plot);
    chartPanel.widgets().set(4, hist_ET);
    chartPanel.widgets().set(5, legendPanel_ndvi); 
    chartPanel.widgets().set(5, legendPanel_et);
 
    // add the output to the map
    //Map.addLayer(ssebi_out.select("ET_daily"), "Daily ET (mm)")
    Map.addLayer(ssebi_out.select('NDVI'),visParamNDVI, "NDVI");
    Map.addLayer(ssebi_out.select('ET_daily'),visParamET, "Daily ET (mm)");
    clearGeometry();
    
}

// var runModel2 = function(){
//   Map.layers().remove("NDVI")
//   Map.layers().remove("Daily ET (mm)")
//   ui.root.remove(controlPanel);
//   ui.root.remove(chartPanel);
//   ui.root.insert(0, controlPanel);
//   ui.root.insert(2, chartPanel)
//   runModel()
// }



// 

var drawingTools = Map.drawingTools();
drawingTools.setShown(false);
function clearGeometry() {
  var layers = drawingTools.layers();
  layers.get(0).geometries().remove(layers.get(0).geometries().get(0));
  Map.setOptions('HYBRID');
}
clearGeometry();

function drawRectangle() {
  clearGeometry();
  drawingTools.setShape('rectangle');
  drawingTools.draw();
}
  
function drawPolygon() {
  clearGeometry();
  drawingTools.setShape('polygon');
  drawingTools.draw();
  }
  


// produce a panel where later the charts will be put
var chartPanel = ui.Panel({
    layout : ui.Panel.Layout.flow("vertical"),
    style:
        {width: '380px', position: 'bottom-right', shown: false}
  })
ui.root.insert(1, chartPanel)
//Map.add(chartPanel);


var RUNbutton = ui.Button({
  label: 'RUN', 
  style: {stretch: 'horizontal', color: '#E48716'}, 
  onClick: function() {
    runModel();
    // Esegui il modello direttamente
  }
});

// drawingTools.onDraw(ui.util.debounce(runModel, 500));
// drawingTools.onEdit(ui.util.debounce(runModel2, 500));

// Panels for the model parameters

var collectionSelector = ui.Select({
  items: [
    { label: 'Landsat 8', value: "LANDSAT/LC08/C02/T1_L2" },
    { label: 'Landsat 9', value: "LANDSAT/LC09/C02/T1_L2" }
  ],
  value: "LANDSAT/LC08/C02/T1_L2",  // Valore di default
  style: {stretch: 'horizontal'}    // Applica lo stretching orizzontale
});

var startDateTexbox = ui.Textbox({
    placeholder: 'startDate (e.g. 2000-01-01)',
    value: '2022-07-01',
    style: { width: '155px', fontSize : '10px', height : '30px' }
});


var RadLimSlider = ui.Slider({
    min: 0.0, max: 0.5, value: 0.11, step: 0.01,
    style: { width: '165px', backgroundColor: bckcolor, color: "blue",
    height : '30px'}
});

var MoistLimSlider = ui.Slider({
    min: 0.0, max: 0.5, value: 0.11, step: 0.01,
    style: { width: '165px', backgroundColor: bckcolor, color: "blue", height : '30px'}
});

var clear=ui.Button({label: 'Clear map' , style: {stretch: 'horizontal', color: '#344C11'}, onClick: function(){
  clearGeometry();
  Map.clear();
  drawingTools.setShown(false);
  Map.setOptions('HYBRID');
  }
});


// build the control panel  
var controlPanel = ui.Panel({
    widgets: [
      ui.Label('This app allows you to estimate and download daily evapotranspiration from Landsat images.', {fontWeight: 'bold', color: '#244072'}),
      ui.Label('1. Select a satellite. Landsat 8 images are available from March 2013, while Landsat 9 images are available from December 2021.', {margin: ' 8px 8px 0px'}),
      ui.Label('2. Select a date. The first available image after the specified date will be displayed.', {margin: '0px 8px'}),
      ui.Label('3. Draw your area of interest by dragging with the mouse.', {margin: '0px 8px'}),
      ui.Label('4. Click on "RUN".', {margin: '0px 8px'}),
      ui.Label('5. Download your image!', {margin: '0px 8px'}),
      ui.Label('If you make a mistake drawing the Area Of Interest, or if you want to start over, click the "Clear map" button and set new parameters.'),
      ui.Label('Define the parameters', {fontWeight: 'bold', margin: 'auto', padding: '10px 0 0 0', color: '#244072'}),
      ui.Label('Select the sensor:'),
      collectionSelector,
      ui.Label('Select the start date:'),
      startDateTexbox,
      ui.Label('Exclude values lower than this albedo threshold for radiation limited line estimation:'),
      RadLimSlider,
      ui.Label('Exclude values lower than this albedo threshold for moisture limited line estimation:'),
      MoistLimSlider,
      ui.Label('Select a drawing mode.'),
      ui.Button({
        label: 'Rectangle',
        onClick: drawRectangle,
        style: {stretch: 'horizontal', height : '40px'}
      }),
      ui.Button({
        label: 'Polygon',
        onClick: drawPolygon,
        style: {stretch: 'horizontal',
          height : '40px'
        }
      }),
      ui.Label('Calculate evapotranspiration', {fontWeight: 'bold', margin: 'auto', padding: '10px 0 0 0', color: '#244072'}),
      RUNbutton,
      clear
    ],
    style: {position: 'bottom-left', width : '380px'},
    layout: ui.Panel.Layout.flow("vertical", true),
  });

// Add the panel to the map
ui.root.insert(0, controlPanel)

  

  
