cd ~/mpgPostdoc/projects/bareEarth/data/canopy/
gdalbuildvrt can.vrt *.tif
gdalwarp can.vrt can.tif
rm can.vrt
gdal_fillnodata.py can.tif canFill.tif
rm can.tif



