# CubeLineMoment.py:

Script to derive Moment0, Moment1, and Moment2 from a set of
input-defined spectral lines in an image cube.  Currently simply
calculates moments over a defined HWZI for each line in band. 

To run in ipython use:

```
run CubeLineMoment.py yaml_scripts/CubeLineMomentInput.yaml
```


## YAML File Input Parameters:

 - `cube` [string]: Input FITS cube to be processed.  Spectral axis can be
   frequency or velocity.
   Example: FITS/NGC253-H2COJ32K02-Feather-line.fits

 - `cuberegion` [string]: ds9 region file used to spatial mask for input FITS
   `cube` emission region.
   Example: regions_folder/NGC253BoxBand6H2COJ32K02.reg

 - `cutoutcube` [string]: Input FITS cube which contains "tracer"
   transition which is strong and representative of dense gas emission
   region traced by other molecules/transitions in `cube`.  Note that this
   `cube` can be any image cube, as long as the PPV range overlaps with `cuberegion`.
   Example: FITS/NGC253-H213COJ32K1-Feather-line.fits

 - `cutoutcuberegion` [string]: ds9 region file used to spatial
   mask input FITS "tracer" transition (`spatialmaskcube`).
   Example: regions_folder/NGC253BoxBand6H2COJ32K02.reg

 - `vz` [float:km/s]: Target central velocity.  In order to maximize the
   effectiveness of the spectral lines extracted from your image cube,
   set `vz` to a value near the median radial velocity of your target.
   Example: 258.8

 - `target` [string]: Target name.  Used as first part of output moment file names.
   Example: NGC253

 - `brightest_line_name` [string]: Brightest line name.
   Example: 13CO_21

 - `brightest_line_frequency` [float:MHz]: Frequency of the bright
   "tracer" transition in `spatialmaskcube`.
   Example: 220.398700

 - `width_line_frequency` [float:MHz]: Frequency of the "representative"
   transition in cube.
   Example: 218.222192

 - `velocity_half_range` [float:km/s]: Estimated half-width at zero
   intensity for the entire velocity extent of the "representative"
   transition in `cube`.  Note that for a galaxy this would be half of
   the total velocity range for the chosen transition.
   Example: 400

 - `noisemapbright_baseline` [list of lists:channels]: Baseline channel segments
   which are considered line-free in `spatialmaskcube`.  Used to
   determine RMS spectral noise in `spatialmaskcube`.
   Example: [[40,60],[100,116],[150,180]]
   
 - `noisemap_baseline` [list of lists:channels]: Baseline channel segments
   which are considered line-free in `cube`.  Used to determine RMS
   spectral noise in `cube`.
   Example: [[20,35],[60,95],[360,370]]

 - `my_line_list` [list:MHz]: List of spectral line frequencies to be
   extracted from `cube`.
   Example: 217.289800, 217.299162, 217.467150, 217.517110, 217.802057, 217.88639, 217.943821, 218.15897, 218.222192, 218.324711, 218.440050, 218.475632, 218.760071, 218.85439, 218.9033555, 218.981019

 - `my_line_widths` [list:km/s]: List of estimated half-width
   zero-intensities for transitions in `my_line_list`.
   Example: 50.0, 50.0, 60.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 50.0, 40.0, 40.0

 - `my_line_names` [list:string]: List of transition names in `my_line_list`.  Used as second part of output moment file names.
   Example: 13CNF122, CH3OH67, 13CNF132, CH3OCHO88, CH3OCHO4847, CH3OH2020, CH3OCHO4546, CH3OCHO??, H2COJ32K0, HC3N2423v0, CH3OH43, H2COJ32K221, H2COJ32K210, HC3N2423v6, OCS1817, HNCO109

 - `signal_mask_limit` [float or None]: Multiplier for noise-based signal
   masking.  Signal less than `signal_mask_limit` times RMS noise is
   masked. Unlike `spatial_mask_limit`, this threshold is applied on a per-voxel basis.
   If this is set to ``None``, no signal masking will be applied.  Use `None` if `min/max_gauss_threshold` masking used.
   Example: 3

 - `spatial_mask_limit` [float or None]: Multiplier for noise-based spatial
   masking.  Signal less than `spatial_mask_limit` times RMS noise is
   masked. Use `None` if `min/max_gauss_threshold` masking used.
   Example: 3

 - `sample_pixel` [string]: File name for ds9 regions file that contains
   the sample pixel positions.  Regions file entries must be of type "point"
   (i.e. point(11.88341,-25.29118) # text={GMC1})
   Example: LeroyNGC253GMCPoint.reg

 - `mask_negatives` [float or bool]: Mask out negatives below N-sigma negative.  Set to `False` if using `min/max_gauss_threshold`.

 - `use_default_width` [bool]: If the width cannot be determined (moment2 is negative, for example), use the `my_line_widths` estimate in place of any pixels with NaN widths.  If using `min/max_gauss_threshold` masking set to `False`, otherwise, set to `True`.

 - `apply_width_mask` [bool]: Should width masking be applied at all?  Turning this off can save some computational time.

 - `min_width` [float:km/s]: The minimum velocity width to allow in the velocity map.  It's a good idea to set this at least to the channel width.

 - `min_gauss_threshold` [float]: Mininum fractional level of `min_width` gaussian to include in noise calculation.  A good value to use is `0.10`.

 - `max_gauss_threshold` [float]: Maximum fractional level of `min_width` gaussian to include in noise calculation.  A good value to use is `0.90`.

 - `use_peak_for_velcut` [bool]: Use the peak velocity to perform the +/- dV velocity cut?  Defaults to `False`, in which case the centroid is used instead.  The centroid is likely more robust, but there are some cases where you might prefer the peak.

 - `dilation_iterations` and `erosion_iterations` [int or None]: Number of itertaions of dilation and erosion to apply to the mask.  This is applied to the width mask and the signal mask with the same number of iterations for both.  Dilation is applied before erosion.  Good values are 4 and 2, respectively.


## Masking Used in CubeLineMoment:

* [optional] Use ds9 regions to select spatial regions to process

* Create a cutout cube `cutoutcube` based on a bright line.
  - [optional] Select only positive values (set by `mask_negatives` parameter)
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

## Do Not Use Moment2

The second moment calculated using this package is biased.  Do not use it.

## Warning Messages

Note that several of these warning messages are due to the use of `NaN` values as blanking values in spectral cubes.  All of these warnings are for information only and can safely be ignored.


```
WARNING: StokesWarning: Cube is a Stokes cube, returning spectral cube for I component [spectral_cube.io.core]
```
Explanation: Cube contains a fourth `Stokes` axis.  Information only.  No action required.


```
WARNING: PossiblySlowWarning: This function (<function BaseSpectralCube.std at 0x19e9053a0>) requires loading the entire cube into memory and may therefore be slow. [spectral_cube.utils]
```
Explanation: `CubeLineMoment` does some calculation on entire cubes.  Information only.  No action required.


```
/Users/jmangum/anaconda3/envs/python39/lib/python3.9/site-packages/numpy/lib/nanfunctions.py:1878: RuntimeWarning: Degrees of freedom <= 0 for slice.
  var = nanvar(a, axis=axis, dtype=dtype, out=out, ddof=ddof,
```
Explanation: There is at least one spectrum that consists of only 1 pixel, so the standard deviation can't be computed. i.e., `noisemask.with_mask(mask[:,None,None]).include().sum(axis=0)` will have at least one pixel with value 1.  Information only.  No action required.


```
/Users/jmangum/anaconda3/envs/python39/lib/python3.9/site-packages/spectral_cube/spectral_cube.py:441: RuntimeWarning: All-NaN slice encountered
```
Explanation: This warning often results from the calculation of the maximum value along the spectral axis toward each pixel in `cutoutcube`.  Since
`cutoutcube` can have blanked (`NaN`) values, there is often at least one position where all spectral values are blanked (`NaN`).  Information only.  No action required.


```
/Users/jmangum/Python/mangum_galaxies-master/CubeLineMoment.py:435: RuntimeWarning: divide by zero encountered in divide
```
Explanation: This warning results from the fact that the denominator in a divide uses a cube with `NaN` values.  Since it is common for a cube to use `NaN` as a blanking value, this warning is common.  Information only.  No action required.


# Comparison Between Masked-Moment and CubeLineMoment Moment0 Calculation:

A detailed comparison using some measurements of HCN 3-2 and 4-3 in NGC253 between a masked-moment algorithm and the algorithm used in CubeLineMoment is described [here](./Moment_Calculation_Comparison.pdf)

# GaussfitGalaxies.py:

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
