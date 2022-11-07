from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Button, Slider, TextInput, Div, Paragraph
from bokeh.layouts import column, row
# from bokeh.io import curdoc
# from bokeh.transform import factor_cmap

from astropy.io import ascii
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, RegularGridInterpolator

# to rename the model, given Gaia phot columns
# for PARSEC models
magRenamer = {
	'G':'phot_g_mean_mag',
	'G_BP' :'phot_bp_mean_mag',
	'G_RP':'phot_rp_mean_mag',
	'g_ps':'g_mean_psf_mag',
	'r_ps':'r_mean_psf_mag',
	'i_ps':'i_mean_psf_mag',
	'z_ps':'z_mean_psf_mag',
	'y_ps':'y_mean_psf_mag',
	'sigG':'phot_g_mean_mag_error',
	'sigG_BP' :'phot_bp_mean_mag_error',
	'sigG_RP':'phot_rp_mean_mag_error',
	'sigg_ps':'g_mean_psf_mag_error',
	'sigr_ps':'r_mean_psf_mag_error',
	'sigi_ps':'i_mean_psf_mag_error',
	'sigz_ps':'z_mean_psf_mag_error',
	'sigy_ps':'y_mean_psf_mag_error',
    'J_2M':'j_m',
	'H_2M':'h_m',
	'Ks_2M':'ks_m',
	'sigJ_2M':'j_msigcom',
	'sigH_2M':'h_msigcom',
	'sigKs_2M':'ks_msigcom'
}

# Redenning coefficients
# from BASE-9 Filters.cpp
absCoeffs0 = {
	"G":      0.86105,
	"G_BP":   1.07185,
	"G_RP":   0.65069,
	"g_ps":   1.16529,
	"r_ps":   0.86813,
	"i_ps":   0.67659,
	"z_ps":   0.51743,
	"y_ps":   0.43092,
	"J_2M":   0.29434,
	"H_2M":   0.18128,
	"Ks_2M":  0.11838,
}
absCoeffs = {}
for key, value in absCoeffs0.items():
	absCoeffs[magRenamer[key]] = value


def getModelGrid(isochroneFile):
	# get the model grid space [age, FeH] from the model file

	FeH = -99.
	age = -99.
	i = 0
	with open(isochroneFile, 'r') as f:
		for line in f:
			if line.startswith('%s'):
				x = line.replace('=',' ').split()
				FeH = float(x[2])
			if line.startswith('%a'):
				x = line.replace('=',' ').split()
				age = float(x[2])

			if (FeH != -99 and age != -99):
				if (line.startswith('%s') or line.startswith('%a')):
					if (i == 0):
						grid = np.array([age, FeH])
					else:
						grid = np.vstack((grid, [age, FeH]))
					i += 1

	return grid


def interpolateModel(age, FeH, isochroneFile, mag = 'phot_g_mean_mag', color1 = 'phot_bp_mean_mag', color2 = 'phot_rp_mean_mag'):
	# perform a linear interpolation in the model grid between the closest points



	# get the model grid
	grid = getModelGrid(isochroneFile)

	# check that the age and Fe/H values are within the grid 
	ageGrid = np.sort(np.unique(grid[:,0]))
	FeHGrid = np.sort(np.unique(grid[:,1]))
	maxAge = np.max(ageGrid)
	minAge = np.min(ageGrid)
	maxFeH = np.max(FeHGrid)
	minFeH = np.min(FeHGrid)

	try:
		# find the 4 nearest age and Fe/H values
		iAge0 = np.where(ageGrid < age)[0][-1]
		age0 = ageGrid[iAge0]
		age1 = age0
		if (iAge0 + 1 < len(ageGrid)):
			age1 = ageGrid[iAge0 + 1]

		iFeH0 = np.where(FeHGrid < FeH)[0][-1]
		FeH0 = FeHGrid[iFeH0]
		FeH1 = FeH0
		if (iFeH0 + 1< len(FeHGrid)):
			FeH1 = FeHGrid[iFeH0 + 1]

		# read in those parts of the isochrone file
		inAge = False
		inFeH = False
		arrS = []
		columns = ''
		testFeH = '-99'
		testAge = '-99'
		with open(isochroneFile, 'r') as f:
			for line in f:
				if (inAge and inFeH and not line.startswith('%') and not line.startswith('#')):
					key = '[' + testAge + ',' + testFeH + ']'
					x = line.strip() + ' ' + testAge + ' ' + testFeH
					arrS.append(x.split())
				if line.startswith('%s'):
					inFeH = False
					x = line.replace('=',' ').split()
					testFeH = x[2]
					if (float(testFeH) == FeH0 or float(testFeH) == FeH1):
						inFeH = True
				if line.startswith('%a'):
					inAge = False
					x = line.replace('=',' ').split()
					testAge = x[2]
					if (float(testAge) == age0 or float(testAge) == age1):
						inAge = True
				if (line.startswith('# EEP')):
					x = line.replace('# ','') + ' logAge' + ' FeH'
					columns = x.split()

		# convert this to a pandas dataframe 
		df = pd.DataFrame(arrS, columns = columns, dtype = float)
		df.rename(columns = magRenamer, inplace = True)

		# take only the columns that we need
		# We cannot interpolate on mass since that is not unique and not strictly ascedning
		# Ideally we would interpolate on initial mass, but that is not included in the BASE-9 model files
		# I will interpolate on EEP, which I *think* is a number that is unique for each star 
		df = df[np.unique(['EEP', mag, color1, color2, 'logAge','FeH'])]
		ages = df['logAge'].unique()
		FeHs = df['FeH'].unique()
		EEPs = df['EEP'].unique()
		
		# initialize the output dataframe
		results = pd.DataFrame()

		# create an array to interpolate on
		# https://stackoverflow.com/questions/30056577/correct-usage-of-scipy-interpolate-regulargridinterpolator
		pts = (ages, FeHs, EEPs)
		for arr in np.unique([mag, color1, color2]):
			val_size = list(map(lambda q: q.shape[0], pts))
			vals = np.zeros(val_size)
			for i, a in enumerate(ages):
				for j, f in enumerate(FeHs):
					df0 = df.loc[(df['logAge'] == a) & (df['FeH'] == f)]
					interp = interp1d(df0['EEP'], df0[arr], bounds_error = False)
					vals[i,j,:] = interp(EEPs)

			interpolator = RegularGridInterpolator(pts, vals)
			results[arr] = interpolator((age, FeH, EEPs))


		return results

	except:
		print(f"!!! ERROR: could not interpolate isochrone.  Values are likely outside of model grid !!!\nage : {age} [{minAge},{maxAge}]\nFeH {FeH} [{minFeH}, {maxFeH}]")
		return None


def createInteractiveIsochrone(photfile, isochroneFile, initialGuess = [4, 0, 0, 0], mag = 'phot_g_mean_mag', color1 = 'phot_bp_mean_mag', color2 = 'phot_rp_mean_mag', xrng = [0.5,2], yrng = [20,10], membershipMin = 0.01):
	'''
	To run this in a Jupyter notebook (intended purpose):
	--------------------
	layout = createInteractiveIsochrone('filename.phot', 'isochrone.model', [log10(age), FeH, mM, Av])
	def bkapp(doc):
		doc.add_root(layout)
	show(bkapp)
	'''

	###########################


	# create the initial figure
	TOOLS = "box_zoom, reset, lasso_select, box_select"
	p = figure(title = "",
		tools = TOOLS, width = 500, height = 700,
		x_range = xrng, y_range = yrng)

	# read in the photometry (this will be replaced by self.data in BASE9_utile)
	data = ascii.read(photfile)

	mask = (data['membership'] > membershipMin) 
	membershipOrg = data['membership'].data.copy() # in case I need to reset

	# add an index column so that I can map back to the original data
	data['index'] = np.arange(0,len(data))

	# get the isochrone at the desired age and metallicity
	iso = interpolateModel(initialGuess[0], initialGuess[1], isochroneFile, mag = mag, color1 = color1, color2 = color2)

	# convert to observed mags given distance modulus and Av for these filters
	# taken from BASE-9 Star.cpp
	# abs is input Av
	# distance is input distance modulus
	# combinedMags[f] += distance;
	# combinedMags[f] += (evoModels.absCoeffs[f] - 1.0) * clust.abs;
	def offsetMag(magVal, magCol, mM, Av):
		return magVal + mM + (absCoeffs[magCol] - 1.)*Av

	def getObsIsochrone(iso, mM, Av):
		color1Obs = offsetMag(iso[color1], color1, mM, Av)
		color2Obs = offsetMag(iso[color2], color2, mM, Av)
		magObs = offsetMag(iso[mag], mag, mM, Av)
		return {'x': color1Obs - color2Obs, 'y': magObs, color1: iso[color1], color2: iso[color2], mag: iso[mag]}



	# define the input for Bokeh
	sourcePhot = ColumnDataSource(data = dict(x = data[mask][color1] - data[mask][color2], y = data[mask][mag], index = data[mask]['index']))
	sourceCluster = ColumnDataSource(data = dict(logAge = [initialGuess[0]], FeH = [initialGuess[1]], mM = [initialGuess[2]], Av = [initialGuess[3]]))
	sourceIso = ColumnDataSource(getObsIsochrone(iso, sourceCluster.data['mM'][0], sourceCluster.data['Av'][0]))

	# add the photometry and isochrone to the plot
	photRenderer = p.scatter(source = sourcePhot, x = 'x', y = 'y', alpha = 0.5, size = 3, marker = 'circle', color = 'black')
	isoRenderer =  p.scatter(source = sourceIso, x = 'x', y = 'y', color = 'red', size = 5)


	###########################
	# widgets
	
	# text boxes to define the age and FeH value for the isochrone
	# adding a bit of extra code to make the label and input side-by-side
	ageInput = TextInput(value = str(initialGuess[0]), title = '')
	ageInputTitle = Paragraph(text = 'log10(age):', align = 'center', width = 60)
	ageInputLayout = row([ageInputTitle, ageInput])
	FeHInput = TextInput(value = str(initialGuess[1]), title = '')
	FeHInputTitle = Paragraph(text = '[Fe/H]:', align = 'center', width = 60)
	FeHInputLayout = row([FeHInputTitle, FeHInput])

	# botton to update the isochrone
	updateIsochroneButton = Button(label = "Update isochrone",  button_type = "success")
	def updateIsochroneCallback(event):
		iso = interpolateModel(float(ageInput.value), float(FeHInput.value), isochroneFile, mag = mag, color1 = color1, color2 = color2)
		if (iso is not None):
			sourceCluster.data['logAge'] = [float(ageInput.value)]
			sourceCluster.data['FeH'] = [float(FeHInput.value)]
			sourceIso.data = getObsIsochrone(iso, sourceCluster.data['mM'][0], sourceCluster.data['Av'][0])
		else:
			ageInput.value = str(sourceCluster.data['logAge'][0])
			FeHInput.value = str(sourceCluster.data['FeH'][0])
	updateIsochroneButton.on_click(updateIsochroneCallback)


	# add sliders to move the isochrone in mM and reddening
	mMSlider = Slider(start = 0, end = 20, value = initialGuess[2], step = 0.01, format = '0.00', title = "Distance Modulus")
	def mMSliderCallback(attr, old, new):
		sourceCluster.data['mM'] = [mMSlider.value]
		iso = sourceIso.data
		sourceIso.data = getObsIsochrone(iso, sourceCluster.data['mM'][0], sourceCluster.data['Av'][0])
	mMSlider.on_change("value", mMSliderCallback)

	AvSlider = Slider(start = 0, end = 3, value = initialGuess[3], step = 0.001, format = '0.000', title = "Av")
	def AvSliderCallback(attr, old, new):
		sourceCluster.data['Av'] = [AvSlider.value]
		iso = sourceIso.data
		sourceIso.data = getObsIsochrone(iso, sourceCluster.data['mM'][0], sourceCluster.data['Av'][0])
	AvSlider.on_change("value", AvSliderCallback)

	# add a reset button
	resetButton = Button(label = "Reset",  button_type = "danger", )

	def resetCallback(event):
		data['membership'] = membershipOrg
		mask = (data['membership'] > membershipMin)
		sourcePhot.data = dict(x = data[mask][color1] - data[mask][color2], y = data[mask][mag], index = data[mask]['index'])

	resetButton.on_click(resetCallback)

	# add a button to write the files
	writeButton = Button(label = "Write .phot and .yaml files",  button_type = "success")

	def writeCallback(event):
		# output updated phot and yaml filesfile
		# self.generatePhotFile()
		# update the yaml variances
		# self.yamlInputDict['Fe_H'] = [self.yamlInputDict['Fe_H'][0], .., ..]
		# self.yamlInputDict['Av'] = [self.yamlInputDict['Av'][0], .., ..]
		# self.yamlInputDict['distMod'] = [self.yamlInputDict['distMod'][0], .., ..]
		# self.yamlInputDict['logAge'] = [self.yamlInputDict['logAge'][0], .., np.inf]
		# self.generateYamlFile()
		# print('Files saved : ', self.photOutputFileName, self.yamlOutputFileName) 
		outfile = photfile + '_new.ecsv'
		ascii.write(data, outfile, overwrite=True) 
		print('File saved : ', outfile) 

	writeButton.on_click(writeCallback)


	# add a button to delete selected points
	deleteButton = Button(label = "Delete selected points",  button_type = "warning")

	def deleteCallback(event):
		# set the membership to -1, redefine the mask, and remove them from the columnDataSource
		if (len(sourcePhot.selected.indices) > 0):
			indices = sourcePhot.data['index'][sourcePhot.selected.indices]
			data['membership'][indices] = -1
			mask = (data['membership'] > membershipMin)
			sourcePhot.data = dict(x = data[mask][color1] - data[mask][color2], y = data[mask][mag], index = data[mask]['index'])
			# reset
			sourcePhot.selected.indices = []

	deleteButton.on_click(deleteCallback)

	###########################
	# layout
	# plot on the left, buttons on the right


	buttons = column(
		Div(text='<div style="height: 15px;"></div>'),
		ageInputLayout,
		FeHInputLayout,
		updateIsochroneButton,
		mMSlider,
		AvSlider,
		deleteButton,
		writeButton,
		Div(text='<div style="height: 50px;"></div>'),
		resetButton,
	)
	title = 	Div(text='<div style="font-size:20px; font-weight:bold">Interactive CMD</div>')
	instructions = 	Div(text='<ul style="font-size:14px">\
		<li>To delete points: select points with lasso or box select tool, and click the "Delete" button.</li>\
		<li>Click the "Reset" button undo all delete actions. </li>\
		<li>To change the isochrone, enter the age and FeH in the appropriate boxes, and click the "Update Isochrone" button.</li>\
		<li>Move the isochrone to change the distance and reddening to better match the data.</li>\
		<li>When finished, click the "Write files" button output the results.</li>\
	</ul>')

	layout = column(title, instructions, row(p, buttons))

	return(layout)



