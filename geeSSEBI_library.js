var bckcolor = "#FFFFFF"
var ERA5H = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY");
var ERA5D = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR");
var scale = 100;
//var visParamET = {bands : ['ET_daily'], max : 5.5, min: 0, opacity :1, palette: ['0af4ff', '39ff51', '98ff4b', 'ffe94d', 'ff7f23', 'ff1004']};
//var visParamNDVI = {bands : ['NDVI'], max : 1, min: 0, opacity :1, palette: ['#DEB887', '#BDB76B', '#8FBC8F', '#556B2F', '#006400']};
var visParamET = {"opacity":1,"bands":["ET_daily"],"min":0,"max":5.5,"palette":['0af4ff', '39ff51', '98ff4b', 'ffe94d', 'ff7f23', 'ff1004']};
var visParamNDVI = {bands : ['NDVI'], max : 1, min: 0, opacity :1, palette: ['#DEB887', '#BDB76B', '#8FBC8F', '#556B2F', '#006400']};
var visParamGeo = {fillColor : '00000000'};



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
  var NDWI = image.select("SR_B3").subtract(image.select("SR_B5"))
  .divide(image.select("SR_B3").add(image.select("SR_B5"))).rename("NDWI")
  return image.addBands(NDWI)
}

var maskWater = function(image){
  return image.updateMask(image.select("NDWI").lt(0.2))
}
    
var getDateProperty = function(feature) {
  return feature.get('date');
};

// Preprocess inputs function

var preprocessInputsL8 = function (startDate,
    endDate,
    geometry) {
    // declare constant
    var B_CONST = ee.Number(5.67).multiply(ee.Number(10).pow(-8));
    //import L8
    var L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filterBounds(geometry)
    var ERA5H = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY");
    var ERA5D = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR");
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
      var NDWI = image.select("SR_B3").subtract(image.select("SR_B5"))
      .divide(image.select("SR_B3").add(image.select("SR_B5"))).rename("NDWI")
      return image.addBands(NDWI)
    }
    
    var maskWater = function(image){
      return image.updateMask(image.select("NDWI").lt(0.2))
    }
        
    var getDateProperty = function(feature) {
      return feature.get('date');
    };
    
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


    return (result)

}

var preprocessInputsL9 = function (startDate,
    endDate,
    geometry) {
    // declare constant
    var B_CONST = ee.Number(5.67).multiply(ee.Number(10).pow(-8));
    //import L8
    var L9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").filterBounds(geometry)
    var ERA5H = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY");
    var ERA5D = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR");
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
      var NDWI = image.select("SR_B3").subtract(image.select("SR_B5"))
      .divide(image.select("SR_B3").add(image.select("SR_B5"))).rename("NDWI")
      return image.addBands(NDWI)
    }
    
    var maskWater = function(image){
      return image.updateMask(image.select("NDWI").lt(0.0))
    }
        
    var getDateProperty = function(feature) {
      return feature.get('date');
    };
    
    // Delcare the cloud masking function and applying scaling factors to the reflectance values
    

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

    var shortwaveNetI = function (image) {
        var result = image.select("ShortwaveDownI").multiply(ee.Image.constant(1).subtract(image.select("albedo"))).rename("ShortwaveNetI")
        return (image.addBands(result))
    }
    var longwaveNetI = function (image) {
        var pt1 = image.select("LongwaveDownI").multiply(image.select('emiss')).subtract(image.select('longwaveEmitted')).rename("LongwaveNetI")//.multiply(image.select("emiss")).rename("LI")
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
        var short = targetRN.select("surface_solar_radiation_downwards_sum")
            .rename("ShortwaveDownD");
        var longwave = targetRN.select("surface_net_thermal_radiation_sum").rename("LongwaveDownD");
        var result = computeNetD(image.addBands(short.clip(geometry)).addBands(longwave.clip(geometry)));
        return (result)
    }

    // final step where we apply the functions declared above
    var an_image = L9.filterDate(startDate, endDate)
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


    return (result)

}

/*
Function prepareInputs
@param startDate: date in format: "YYYY-MM-DD"
@param endDate: date in format: "YYYY-MM-DD" (shall be greater than start date)
@param geometry: object of class ee.Geometry; a polygon identifying the study area.
*/

var preprocessInputsL5 = function(startDate,
                              endDate,
                              geometry){
  // declare constant
  var B_CONST = ee.Number(5.67).multiply(ee.Number(10).pow(-8));
  // import Landsat 5 image collection and ERA5 Land both hourly and daily
  var L5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterBounds(geometry)
  var ERA5H = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY");
  var ERA5D = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR");
  // Delcare the cloud masking function and applying scaling factors to the reflectance values
    var maskL5sr = function(image) {
    // Bit 0 - Fill
    // Bit 1 - Dilated Cloud
    // Bit 2 - Cirrus
    // Bit 3 - Cloud
    // Bit 4 - Cloud Shadow
    var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
    var saturationMask = image.select('QA_RADSAT').eq(0);
    
    //My own cloud distance mask to be more restrictive
    var cdist = image.select("ST_CDIST").multiply(0.01);
    var cdist_mask = cdist.gt(ee.Image.constant(10));
    
    // Apply the scaling factors to the appropriate bands.
    var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
    var thermalBands = image.select('ST_B.').multiply(0.00341802).add(149.0);
    var constant = ee.Image.constant(0.0001).rename('area');
    // Aggiungere maschera per togliere valori di ST eccessivi, tipo maggiori di 350 K e minori di 273?
    // Replace the original bands with the scaled ones and apply the masks.
    return image.addBands(opticalBands, null, true)
      .addBands(thermalBands, null, true)
      .addBands(constant)
      .updateMask(qaMask)
      .updateMask(saturationMask);
  };
  // compute NDWI for water masking
  var computeNDWI = function(image){
    var NDWI = image.select("SR_B2").subtract(image.select("SR_B5"))
    .divide(image.select("SR_B2").add(image.select("SR_B5"))).rename("NDWI");
    return image.addBands(NDWI);
  }
  
  var maskWater = function(image){
    return image.updateMask(image.select("NDWI").lt(0.0))
  }
  // Declare functions for computing NDVI and other spectral properties
  // They are normalized to be restricted into the range between 0 and 1
  
  function obtainNDVI(image){
    var NDVI = image.select('SR_B4').subtract(image.select('SR_B3'))
    .divide(image.select('SR_B4').add(image.select('SR_B3')))
    .rename('NDVI');
    var  NDVI_ok = NDVI.multiply(NDVI.gt(ee.Image.constant(0))).multiply(NDVI.lt(ee.Image.constant(1)))
    .add(NDVI.gt(ee.Image.constant(1))).rename("NDVI")
    var result = image.addBands(NDVI_ok);
    return(result)
  }


  function obtainFC(image){
    var FC = image.select("NDVI").subtract(ee.Image.constant(0.2))
    .divide(ee.Image.constant(0.6)).pow(ee.Image.constant(2)).rename("FC");
    FC = FC.where(FC.gt(1), ee.Image.constant(1)).where(image.select("NDVI").lt(0.2), ee.Image.constant(0))
    var result = image.addBands(FC);
    return(result)
  }
  
  function obtainALBEDO(image){
    // Landsat 5 albedo: 0.356*B1 + 0.130 B3 + 0.373 * B4 + 0.085*B5 + 0.072 * B7 - 0.0018
    // Liang 2001 https://doi.org/10.1016/S0034-4257(00)00205-4
    var ALBEDO = image.select('SR_B1').multiply(ee.Image.constant(0.356))
    .add( image.select('SR_B3').multiply(ee.Image.constant(0.130)))
    .add( image.select('SR_B4').multiply(ee.Image.constant(0.373)))
    .add( image.select('SR_B5').multiply(ee.Image.constant(0.085)))
    .add( image.select('SR_B7').multiply(ee.Image.constant(0.072)))
    .add(ee.Image.constant(-0.0018))
    .rename('albedo');
    var  ALBEDO_ok = ALBEDO.multiply(ALBEDO.gt(ee.Image.constant(0))).multiply(ALBEDO.lt(ee.Image.constant(1)))
    .add(ALBEDO.gt(ee.Image.constant(1)))
    var result = image.addBands(ALBEDO_ok.rename('albedo'));
    return(result)
  }
  //rewrite and add
function obtainEMISS(image){
  var emiss = ee.Image.constant(0.971).multiply(ee.Image.constant(1).subtract(image.select("FC")))
  .add(ee.Image.constant(0.982).multiply(image.select("FC"))).rename("emiss");
  var result = image.addBands(emiss);
  return(result);
}
  /* ========== Preprocess function ============
  Produces the albedo discrete band, for the computation of the dry and wet regression lines
  Renames bands and selects the bands necessary for running the SEBS model.
  */
  var preprocess = function(image){
    var albedo_disc = image.select('albedo').multiply(ee.Image.constant(100)).round().rename('albedo_disc');
    var result = image.select('albedo')
    .addBands(image.select("NDVI"))
    .addBands(image.select("FC"))
    .addBands(image.select('emiss'))
    .addBands(image.select('ST_B6').rename('LST'))
    .addBands(albedo_disc);
    return(result);
  };
  
  // function for emitted radiation
  var computeEmitted = function(image){
    var emitted = image.select('LST').pow(4);
    emitted = emitted.multiply(image.select('emiss')).multiply(ee.Image.constant(B_CONST)).rename("longwaveEmitted");
    var result = image.addBands(emitted);
    return(result);
  };
  

  var shortwaveNetI = function(image){
    var result = image.select("ShortwaveDownI").multiply(ee.Image.constant(1).subtract(image.select("albedo"))).rename("ShortwaveNetI")
    return(image.addBands(result))
  }
  var longwaveNetI = function(image){
    var pt1 = image.select("LongwaveDownI").multiply(image.select('emiss')).subtract(image.select('longwaveEmitted')).rename("LongwaveNetI")//.multiply(image.select("emiss")).rename("LI")
    //var pt2 = image.select("emiss").multiply(image.select("LST").pow(ee.Image.constant(4))).multiply(ee.Image.constant(B_CONST)).rename("LO")
    //var mean_emitted = pt2.reduceRegion(ee.Reducer.mean(),geometry, 100, "EPSG:32632")
    //print("Longwave emitted")
    //print(mean_emitted.get("LO"))
    //var result = pt1.select("LI").subtract(pt2.select("LO")).rename("LongwaveNetI")
    return(image.addBands(pt1))
  }
  var shortwaveNetD = function(image){
    var result = image.select("ShortwaveDownD").multiply(ee.Image.constant(1).subtract(image.select("albedo"))).rename("ShortwaveNetD")
    return(image.addBands(result))
  }
  var longwaveNetD = function(image){
    var pt1 = image.select("LongwaveDownD").rename("LongwaveNetD")//.multiply(image.select("emiss"))
    //.rename("LI")
    //var pt2 = image.select("emiss").multiply(image.select("LST").pow(ee.Image.constant(4))).multiply(ee.Image.constant(B_CONST)).rename("LO")
    //var result = pt1.select("LI").subtract(pt2.select("LO")).rename("LongwaveNetD")
    return(image.addBands(pt1))
  }
  
  var computeNetI = function(image){
    var result = shortwaveNetI(image).select("ShortwaveNetI").add(longwaveNetI(image).select("LongwaveNetI"))
    .rename("RnI")
    return(image.addBands(result))
  }
    
  var computeNetD = function(image){
    var result = shortwaveNetD(image).select("ShortwaveNetD").multiply(ee.Image.constant(1).subtract(image.select('albedo'))).add(longwaveNetD(image).select("LongwaveNetD"))
        .rename("RnD")
    return(image.addBands(result))
  }
  
  
  var obtainRadiationI = function(image,
                                  startDate,
                                  endDate,
                                  hour){
    // migliorare il filtraggio per prendere solo quella con timestamp esatto uguale all'immagine L8 
    // anche quando si dovesse usare un range piu lungo
    var targetRN = ERA5H.filterDate(image.date(), image.date().advance(1, "day")).filter(ee.Filter.eq("hour", hour)).first();
    //print("hourly ERA5")
    //print(targetRN)
    var short = targetRN.select("surface_solar_radiation_downwards_hourly")
    .divide(ee.Image.constant(3600)).rename("ShortwaveDownI");// dividi per 3600 per ottenere istantanea
    
    var longwave = targetRN.select("surface_thermal_radiation_downwards_hourly").divide(ee.Image.constant(3600)).rename("LongwaveDownI");
    var result = computeNetI(image.addBands(short.clip(geometry)).addBands(longwave.clip(geometry)));
    
    
    return(result)
  }
  var obtainRadiationD = function(image, 
                                  startDate,
                                  endDate){
                                    
    // migliorare il filtraggio per prendere solo quella con timestamp esatto uguale all'immagine L8 
    // anche quando si dovesse usare un range piu lungo
    var targetRN = ERA5D.filterDate(image.date(),image.date().advance(1, "day")).first();
    //print(targetRN)
    var short = targetRN.select("surface_solar_radiation_downwards_sum")
    .rename("ShortwaveDownD");
    var longwave = targetRN.select("surface_net_thermal_radiation_sum").rename("LongwaveDownD");
    var result = computeNetD(image.addBands(short.clip(geometry)).addBands(longwave.clip(geometry)));
    return(result)
  }
  
    // final step where we apply the functions declared above
  var an_image = L5.filterDate(startDate, endDate)
  .filterBounds(geometry).map(maskL5sr)
  .map(computeNDWI).map(maskWater)
  .map(obtainALBEDO)
  .map(obtainNDVI).map(obtainFC).map(obtainEMISS)
  .map(preprocess).map(computeEmitted)
  .first().clip(geometry);
  var time = an_image.date();// get the time
  var hour = time.millis().divide(60 * 60 * 1000 * 24) // convert time to millis
    .subtract(time.millis().divide(60 * 60 * 1000 * 24).round()).multiply(24).round() //subtract the day to get the hour

  var result = obtainRadiationD(obtainRadiationI(an_image, 
                                startDate, 
                                endDate,
                                hour),
                                startDate, 
                                endDate)
                                

  return(result)
                                
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
    var CRS = image.projection().crs()
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
        //print(hist.get('albedo_disc'))
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
        var result = image.addBands(EF.rename('EF'))
            .addBands(ee.Image.constant(ee.Number(yIntLR)).rename('yIntLR'))
            .addBands(ee.Image.constant(ee.Number(slopeLR)).rename('slopeLR'))
            .addBands(ee.Image.constant(ee.Number(yIntHR)).rename('yIntHR'))
            .addBands(ee.Image.constant(ee.Number(slopeHR)).rename('slopeHR'));
        return (result)
    }
    // Finally compute
    var provaET = computeET(addSFCoef(computeEF(image)));
    provaET = provaET.addBands(ee.Image.constant(provaET.date().millis().divide(1000)).rename('date_seconds'))
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
                      title: 'albedo',
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
// function that produces the ET histogram
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
              title: 'frequency',
              titleTextStyle: {italic: false, bold: true}
            }
          })
    return(et_chart)
}


exports.preprocessInputsL5 = preprocessInputsL5;
exports.preprocessInputsL8 = preprocessInputsL8;
exports.preprocessInputsL9 = preprocessInputsL9;
exports.runSSEBI = runSSEBI;
exports.plotSSEBI = plotSSEBI;
exports.etHist = etHist;
exports.visParamET = visParamET;
exports.visParamNDVI = visParamNDVI;

