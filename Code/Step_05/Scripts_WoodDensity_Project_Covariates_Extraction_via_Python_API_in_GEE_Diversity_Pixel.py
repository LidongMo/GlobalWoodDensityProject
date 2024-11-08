import ee
from time import sleep
import multiprocessing
import math
import time
import sys
import pandas as pd
import numpy as np 
from pathlib import Path
from functools import partial
from contextlib import contextmanager
from ctypes import c_int
from multiprocessing import Value, Lock, Process

# python3 Scripts_WoodDensity_Project_Covariates_Extraction_via_Python_API_in_GEE_Diversity_Pixel.py  >> output.txt

ee.Initialize()

# Composite image to sample, you can define you own composite and put it at google earth engine Asset!
compositeToUseRaw = ee.Image("users/leonidmoore/ForestBiomass/20200915_Forest_Biomass_Predictors_Image")
# get the projection
stdProj = compositeToUseRaw.projection()
# load the age bands and rename it
forestAgeData = ee.Image("projects/crowtherlab/johan/Besnard_ForestAge")
# we choose the band "forest_age_TC030" for following modeling
forestAge = forestAgeData.select(['forest_age_TC000']).reproject(crs=stdProj).rename('ForestAge')
# load the additional layers and uniform the projection
soilmoisture = ee.Image('users/haozhima95/wld_soil_moisture').reproject(crs=stdProj).rename('SoilMoisture')
aridityIndexVer3 = ee.Image('users/leonidmoore/Aridiy_index_Ver_3_2023').rename('AridityIndex3').reproject(stdProj)
ET_Ver3 = ee.Image('users/leonidmoore/Evapotranspiration_0_Ver_3_2023').rename('ET0').reproject(stdProj)
# add the bands to the composite
compositeToUse = compositeToUseRaw.addBands(soilmoisture).addBands(forestAge).addBands(aridityIndexVer3).addBands(ET_Ver3)

# FeatureCollection to sample. the path is the right path of the feature collection in Google earth engine Asset
# points = ee.FeatureCollection("users/leonidmoore/WoodDensityProject/Full_Shapefiles/WoodDensity_Diversity_Pixel") this is the one without the 
points = ee.FeatureCollection("users/leonidmoore/WoodDensityProject/Full_Shapefiles/WoodDensity_Diversity_SpNr_updated_2023")
print('Number of points to sample:', points.size().getInfo())

def FCDFconv(fc):
		features = fc.getInfo()['features']
		dictarray = []

		for f in features:
				dict = f['properties']
				dictarray.append(dict)

		df = pd.DataFrame(dictarray)

		return df

origCols = list(FCDFconv(points.limit(1)).columns)

# Bands to sample. Default: all bands plus variables present in original FC
BANDS = compositeToUse.bandNames().getInfo()+origCols

counter = multiprocessing.Value(c_int)  # defaults to 0
counter_lock = Lock()

def increment(results):
		with counter_lock:
				counter.value += len(results)

def extract_grid(region, points):
		"""
		Extracts a single point.
		This handles the too-many-requests error by idling the worker with backoff.
		"""
		success = False
		idle = 0

		result = []
		while not success:
				try:
						if int(points.filterBounds(ee.Feature(region).geometry()).size().getInfo()) < 12000:

									values = (compositeToUse.reduceRegions(collection = points.filterBounds(ee.Feature(region).geometry()),
																												reducer = ee.Reducer.first(),
																												scale = compositeToUse.projection().nominalScale().getInfo(),
																												tileScale = 16)
														.toList(50000)
														.getInfo())

									for item in values:
											values = item['properties']
											row = [str(values[key]) for key in BANDS]
											row = ",".join(row)
											result.append(row)

									if len(result) > 0:
										print("Processed %d features" % len(result))

									return result

						else:
								pointsWithRandom = points.filterBounds(ee.Feature(region).geometry()).randomColumn()

								nPoints = pointsWithRandom.size().getInfo()

								# Number of subsets when trying to sample a FC that is too large
								nSubsets = round(int((nPoints/9000)), 0)+1

								# List of bins
								mapList = list(range(1, nSubsets+1))

								# List of breakpoints
								breakPoints = [x / nSubsets for x in [0]+mapList]

								print('FC too large: ', nPoints, ', splitting in ', nSubsets, ' subsets')

								result_sub = []
								# success = False

								for i in mapList:
										result_sub = []

										values = (compositeToUse.reduceRegions(collection = pointsWithRandom.filter(ee.Filter.And(
																																											ee.Filter.gte('random', breakPoints[i - 1]),
																																											ee.Filter.lte('random', breakPoints[i]))),
																																	reducer = ee.Reducer.first(),
																																	scale = compositeToUse.projection().nominalScale().getInfo(),
																																	tileScale = 16)
																									.toList(50000)
																									.getInfo())

										for item in values:
													values = item['properties']
													row = [str(values[key]) for key in BANDS]
													row = ",".join(row)
													result_sub.append(row)

										# print('Items added in subset ', i, ":", len(result_sub))

										result = result + result_sub

										print('Processed ', len(result), ' from', nPoints, 'after ', i, '/', nSubsets, 'subsets')

								return result

				except Exception as e:
							print(e)
							success = False
							idle = (1 if idle > 5 else idle + 1)
							print("idling for %d" % idle)
							sleep(idle)

def extract_and_write_grid(n, grids, points, region):
    region = grids.get(n)
    results = extract_grid(region, points)

    if len(results) > 0:
        print("Processed %d features" % len(results))

        df = pd.DataFrame([item.split(",") for item in results], columns = BANDS)
        df.replace('None', np.nan, inplace = True)
        df.to_csv("CovariatesExtractionFolder/SubExtractedTables/sampled_%d.csv" % n)

def generateGrid(region, size):
	"""Generate a grid covering the region with size*size rectangles"""
	bins = ee.Number(size)
	coords = ee.List(region.coordinates().get(0))
	xs = coords.map(lambda l : ee.List(l).get(0))
	ys = coords.map(lambda l : ee.List(l).get(1))

	xmin = ee.Number(xs.reduce(ee.Reducer.min()))
	xmax = ee.Number(xs.reduce(ee.Reducer.max()))
	ymin = ee.Number(ys.reduce(ee.Reducer.min()))
	ymax = ee.Number(ys.reduce(ee.Reducer.max()))

	dx = xmax.subtract(xmin).divide(bins)
	dy = ymax.subtract(ymin).divide(bins)

	def f1(n):
		def f2(m):
			x1 = xmin.add(dx.multiply(n))
			y1 = ymin.add(dy.multiply(m))
			return ee.Geometry.Rectangle([x1, y1, x1.add(dx), y1.add(dy)], None, False)
		return ee.List.sequence(0, bins.subtract(1)).map(f2)
	grid = ee.List.sequence(0, bins.subtract(1)).map(f1).flatten().flatten()
	return grid

@contextmanager
def poolcontext(*args, **kwargs):
		"""This just makes the multiprocessing easier with a generator."""
		pool = multiprocessing.Pool(*args, **kwargs)
		yield pool
		pool.terminate()

if __name__ == '__main__':
		unboundedGeo = ee.Geometry.Polygon([[[-180, 88], [180, 88], [180, -88], [-180, -88]]], None, False);

		terra_poly = compositeToUse.select(BANDS[0]).toInt().reduceToVectors(
												geometry= unboundedGeo,
												crs= compositeToUse.select(BANDS[0]).projection(),
												scale= compositeToUse.select(BANDS[0]).projection().nominalScale().getInfo(),
												geometryType= 'polygon',
												eightConnected= False,
												labelProperty= 'zone',
												# reducer= ee.Reducer.mean(),
												maxPixels= 1E13
											);

		GRID_WIDTH = 40  # How many grid cells to use.
		grids = generateGrid(unboundedGeo, GRID_WIDTH)
		size = grids.size().getInfo()
		print("Grid size: %d " % size)

		# How many concurrent processors to use.  If you're hitting lots of
		# "Too many aggregation" errors (more than ~10/minute), then make this
		# number smaller.  You should be able to always use at least 20.
		NPROC = 50
		with poolcontext(NPROC) as pool:
			results = pool.map(
				partial(extract_and_write_grid, grids=grids, points=points, region=terra_poly),
				range(0, size))
