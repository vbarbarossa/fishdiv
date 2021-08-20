# on RU server
# cropland = 10 and 20

for(year in c(1992,2015)){
  
  system(
    paste0(
      'gdal_calc.py -A /vol/milkundata/ESA_landcover/TIFF/ESACCI-LC-L4-LCCS-Map-300m-P1Y-',year,'-v2.0.7.tif',
      ' --calc="logical_or(A == 10,A == 20)" --outfile=cropland_',year,'.tif --type=Byte --co="COMPRESS=LZW"'
    )
  )
}

