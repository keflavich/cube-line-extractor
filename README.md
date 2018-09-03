# mangum_galaxies

## CubeLineMoment.py:

Script to derive Moment0, Moment1, and Moment2 from a set of
input-defined spectral lines in an image cube.  Currently simply
calculates moments over a defined HWZI for each line in band. 

To run in ipython use:

run CubeLineMoment.py yaml_scripts/CubeLineMomentInput.yaml


YAML File Input Parameters:

-- cube [string]: Input FITS cube to be processed.  Spectral axis can be
   frequency or velocity.
   Example: FITS/NGC253-H2COJ32K02-Feather-line.fits

-- cuberegion [string]: ds9 region file used to spatial mask for input FITS
   cube emission region.
   Example: regions/NGC253BoxBand6H2COJ32K02.reg

-- cutoutcube [string]: Input FITS cube which contains "tracer"
   transition which is strong and representative of dense gas emission
   region traced by other molecules/transitions in cube.
   Example: FITS/NGC253-H213COJ32K1-Feather-line.fits

-- cutoutcuberegion [string]: ds9 region file used to spatial
   mask input FITS spatialmaskcube.
   Example: regions/NGC253BoxBand6H2COJ32K02.reg

-- vz [float:km/s]: Target central velocity.  In order to maximize the
   effectiveness of the spectral lines extracted from your image cube,
   set vz to a value near the median radial velocity of your target.
   Example: 258.8

-- target [string]: Target name.
   Example: NGC253-H2COJ32K02

-- brightest_line_frequency [float:MHz]: Frequency of the bright
   "tracer" transition in spatialmaskcube.
   Example: 219.560358

-- width_line_frequency [float:MHz]: Frequency of the "representative"
   transition in cube.
   Example: 218.222192

-- velocity_half_range [float:km/s]: Estimated half-width at zero
   intensity for the entire velocity extent of the "representative"
   transition in cube.  Note that for a galaxy this would be half of
   the total velocity range for the chosen transition.
   Example: 80

-- noisemapbright_baseline [list of lists:channels]: Baseline channel segments
   which are considered line-free in spatialmaskcube.  Used to
   determine RMS spectral noise in spatialmaskcube.
   Example: [[40,60],[100,116],[150,180]]
   
-- noisemap_baseline [list of lists:channels]: Baseline channel segments
   which are considered line-free in cube.  Used to determine RMS
   spectral noise in cube.
   Example: [[20,35],[60,95],[360,370]]

-- my_line_list [list:MHz]: List of spectral line frequencies to be
   extracted from cube.
   Example: 217.289800, 217.299162, 217.467150, 217.517110, 217.802057, 217.88639, 217.943821, 218.15897, 218.222192, 218.324711, 218.440050, 218.475632, 218.760071, 218.85439, 218.9033555, 218.981019

-- my_line_widths [list:km/s]: List of estimated half-width
   zero-intensities for transitions in my_line_list.
   Example: 50.0, 50.0, 60.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 50.0, 40.0, 40.0

-- my_line_names [list:string]: List of transition names in my_line_list.
   Example: 13CNF122, CH3OH67, 13CNF132, CH3OCHO88, CH3OCHO4847, CH3OH2020, CH3OCHO4546, CH3OCHO??, H2COJ32K0, HC3N2423v0, CH3OH43, H2COJ32K221, H2COJ32K210, HC3N2423v6, OCS1817, HNCO109

-- signal_mask_limit [float]: Multiplier for noise-based signal
   masking.  Signal less than signal_mask_limit times RMS noise is
   masked. 
   Example: 2

-- spatial_mask_limit [float]: Multiplier for noise-based spatial
   masking.  Signal less than spatial_mask_limit times RMS noise is
   masked. 
   Example: 2



Masking Used in CubeLineMoment:

-- [optional] Using an input ds9 regions file, mask only those spatial
regions containing emission to process


-- Using a limit on signal intensity derived from the noise in
spatialmaskcube (which is defined in the parameter noisemapbright).
Algorithm is any noisemapbright pixel with intensity (in absolute
value) which is greater than spatial_mask_limit * noisemapbright.

-- Define a slab which uses spatialmaskcube (the designated “brightest
line” which defines where we expect to find dense gas) which has been
slabbed to a velocity width of +-linewidth_guess around the brightest
line frequency and spatially masked by spatialmaskcuberegion (another
ds9 file, which is usually the same as cuberegion.

-- Define a subcube which is a velocity slab of cube which has velocity
width peak_velocity.min() - line_width to peak_velocity_max() +
line_width.  The peak_velocity list is derived from the peak
velocities from the designated “brightest line” in brightest_cube.

-- For each transition to be evaluated, define a gaussian mask
(gauss_mask_cube) centered on the centroid (first moment) for the
transition used to define the standard line width
(width_line_frequency). 

-- Define a gaussian threshold for the gaussian which scales as the
peak-signal-to-noise, or 1/peak-signal-to-noise.

-- Apply the threshold to the gaussian mask such that all pixels with
intensity greater than 1/peak-signal-to-noise are passed:
width_mask_cube = gauss_mask_cube > threshold.

-- Define a mask which includes only spectral planes with velocities 
within a linewidth of the peak velocity.

-- Define a signal_mask using the noisemap defined above:
signal_mask = subcube > signal_mask_limit * noisemap.

-- Apply remaining masks: mask, spatial_mask, signal_mask.

# Masking

* [optional] Use ds9 regions to select spatial regions to process
   (This should not be used, since it is not supported in later steps)

* Create a cutout cube `cutoutcube` based on a bright line.
  - [optional] Select only positive values
  - Select a subset of the cube at +/- `velocity_half_range` from 
    the central velocity `vz`
  - Compute peak intensity `max_map`, width `width_map`, and peak velocity
    `peak_velocity` of this cube to use in future steps

* Create a noise map `noisemapbright` based on the bright line cube
  - Select signal-free baseline regions using the `noisemapline_baseline` parameter
  - Compute the standard deviation in the spectral direction of the selected region

* Create another noise map `noisemap` based on the target cube (the process is
  the same as for the bright cube)

* Create a spatial mask based on the peak intensity of the bright line cutout cube `cutoutcube`:
  pixels in the peak map `max_map` of the cutout cube above `signal_mask_limit` *
  `noisemapbright` are included

* Using the bright line maps, make a Gaussian mask cube `gauss_mask_cube` for each target line
 - For each included _spatial_ pixel, produce a Gaussian spectrum using the centroid
 from `peak_velocity`, the peak intensity from `max_map`, and the width from `width_map`
 - Compute the peak signal-to-noise in each pixel by taking `max_map` / `noisemap`
 - Determine a threshold that is 1/`peak_sn`
 - Create a PPV inclusion mask `width_mask_cube` where `gauss_mask_cube` > `threshold`

* [optional] Create a S/N mask where any PPV pixel is greater than
  `signal_mask_limit` * `noisemap` (this is a comparison between a cube and a
  spatial map)

* Create a PPV mask `velocity_range_mask` where the velocity is within
  `line_width` of `peak_velocity`

* Select the data combining the `velocity_range_mask`, the S/N limit, and the
  Gaussian-based `width_mask_cube`


## GaussfitGalaxies.py:

An implementation of gaussfit_catalog (see
https://github.com/radio-astro-tools/gaussfit_catalog) using astropy models to
do gaussian fits to a list of input FITS files using a list of input positions
from a DS9 regions file. 

If you have input regions which are off the image, the script will squawk but not crash.

The input regions are the initial guesses, but they should be very close to the
peak.  Note that gaussfit_catalog is catered to situations where the background
was a significant confusing factor once you got more than 1-2 beam FWHM from
the peak, so it is quite restrictive in how far it will wander beyond the
initial position guess.

The four panels in the output png files are showing:
* Top Left: Data
* Top Right: Fit
* Bottom Left: Residual
* Bottom Right: Data with the fit contoured on top

The .ipac files are in the "ascii.ipac” format from astropy.  In those files,
the non-obvious columns are: 
* chi2_n = chi2/n_pixels, which is close to a reduced chi2
* PA is defined as east-from-north (but double check this!  angle conventions
  are tricky)

Two examples of a similar implementation of gaussfit_catalog are:
https://github.com/keflavich/W51_ALMA_2013.1.00308.S/blob/0789ccbb2fd3bfe801cfb63818ad2696825d076f/analysis/longbaseline/gaussfit_sources.py
https://github.com/keflavich/SgrB2_ALMA_3mm_Mosaic/blob/93fe253f6de499cc91779bdae4f0e22ab806c161/analysis/gaussfit_sources.py
