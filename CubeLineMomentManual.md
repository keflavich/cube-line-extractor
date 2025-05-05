
# CubeLineMoment

Script to derive Moment0, Moment1, and Moment2 from a set of input-defined spectral lines in an image cube. Currently simply calculates moments over a defined HWZI for each line specified.

## N.B.
1. It is important to keep in mind that one needs to set `sample_pixel` in order to get diagnostic plots, otherwise 
none will be produced. 
2. The **choice of sample pixel(s)** is significant, and can impact the quality of the moment maps generated.
It is best to choose one or more of the strongest peaks in the intensity for diagnostic calculation and display (see below).
3. Regarding the choise of `cutoutcube`, if you are extracting moments from relative strong spectral line(s), setting `cutoutcube` equal to `cube` is the best option.  If, on the other hand, your spectral line(s) of interest are weak, choose a `cutoutcube` with a bright spectral line that can be used to guide the extraction of moments for your spectral line(s) of interest.

## Requirements

    aplpy
    pyspeckit
    spectral-cube
    radio-beam
    yaml

To run in ipython use:
```sh
% run CubeLineMoment.py yaml_scripts/CubeLineMomentInput.yaml
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


## Masking Used in CubeLineMoment

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


## Worked Example

In the following we will show how a typical run of **CubeLineMoment** will look.  For this example I am processing an ALMA Band 6 measurement of the galaxy NGC253 from projects 2018.1.00765.S, 2021.1.00105.S, and 2019.1.00113.S.  Once I have edited my **yaml** file appropriately, it looks like the following:
```python
cube: /Users/jmangum/Science/ALCHEMI/ALCHEMIHighResBand67/Imaging/ngc253_h_06_TM1/spw29/uid___A001_X133d_X1579.ngc253_sci.spw29.cube.selfcal.I.fits
cuberegion: CubeLineMoment.reg
cutoutcube: /Users/jmangum/Science/ALCHEMI/ALCHEMIHighResBand67/Imaging/ngc253_h_06_TM1/spw29/uid___A001_X133d_X1579.ngc253_sci.spw29.cube.selfcal.I.fits
cutoutcuberegion: CubeLineMoment.reg
vz: 236
target: NGC253
brightest_line_name: HCN_32
brightest_line_frequency: 265.88618
width_line_frequency: 265.88618
velocity_half_range: 250
noisemapbright_baseline: [[55,60],[150,165]]
noisemap_baseline: [[55,60],[150,165]]
my_line_list: 265.88618
my_line_widths: 120
my_line_names: HCN_32
signal_mask_limit: None
spatial_mask_limit: None
mask_negatives: False
sample_pixel: Region6.reg
use_default_width: True
min_width: 10
min_gauss_threshold: 0.10
max_gauss_threshold: 0.90
dilation_iterations: 4
erosion_iterations: 2
```
...where you can see that I am using the HCN 3-2 transition as a my bright "tracer" and "target" transitions.  The **CubeLineMoment.reg** and **Region6.reg** regions files look like the following...
```python
# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite-1 dash=0 fixed=0 edit=1 move=1 delete=1 source=1 include=1
fk5
box(00:47:33.120,-25:17:17.59,60.0",50.0",0)
```

```python
 Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
point(0:47:33.3312,-25:17:15.756) # text={Region 6}
```

I can run **CubeLineMoment** as follows in ipython:

```python
In [2]: run CubeLineMoment.py yaml_scripts/HCN32.yaml
WARNING: StokesWarning: Cube is a Stokes cube, returning spectral cube for I component [spectral_cube.io.core]
WARNING: StokesWarning: Cube is a Stokes cube, returning spectral cube for I component [spectral_cube.io.core]
WARNING: PossiblySlowWarning: This function (<function BaseSpectralCube.argmax at 0x10b97b060>) requires loading the entire cube into memory and may therefore be slow. [spectral_cube.utils]
WARNING: PossiblySlowWarning: This function (<function BaseSpectralCube.max at 0x10b97aca0>) requires loading the entire cube into memory and may therefore be slow. [spectral_cube.utils]
/Users/jmangum/anaconda3/envs/py313/lib/python3.13/site-packages/spectral_cube/spectral_cube.py:436: RuntimeWarning: All-NaN slice encountered
  out = function(self._get_filled_data(fill=fill,
/Users/jmangum/anaconda3/envs/py313/lib/python3.13/site-packages/numpy/lib/_nanfunctions_impl.py:2019: RuntimeWarning: Degrees of freedom <= 0 for slice.
  var = nanvar(a, axis=axis, dtype=dtype, out=out, ddof=ddof,
noisemapbright peak = 0.0013851320836693048 Jy / beam
/Users/jmangum/anaconda3/envs/py313/lib/python3.13/site-packages/numpy/lib/_nanfunctions_impl.py:2019: RuntimeWarning: Degrees of freedom <= 0 for slice.
  var = nanvar(a, axis=axis, dtype=dtype, out=out, ddof=ddof,
INFO: There are 341655 bad (nan) values in the width map [__main__]

INFO: Line: HCN_32, 265.88618 GHz, 120.0 km / s [__main__]
Peak S/N: 375.266845703125
Minimum S/N: -1.6258751153945923
Highest Threshold: 0.8999999761581421
Lowest Positive Threshold: 0.10000000149011612
Lowest Threshold: 0.10000000149011612
Sample Pixel =  (375, 436, 'Region 6')
SP Threshold: 0.10000000149011612
SP S, N, S/N: 0.015526569448411465 Jy / beam, 0.0007053805748000741 Jy / beam, 22.011619567871094
Number of values above threshold: 14501727 = 1.45e+07
Number of spatial pixels excluded: 133391 out  of 613089
Number of spatial pixels excluded in spatially included region: 2021 out  of 481719
Min, Max value in the mask cube: 0.0,0.9999999999999996
shapes: mask cube=(74, 783, 783)  threshold: (783, 783)
WARNING: AstropyDeprecationWarning: CoordinateHelper.ticks should not be accessed directly and is deprecated [astropy.visualization.wcsaxes.coordinate_helpers]
INFO: Auto-setting vmin to -9.594e-01 [aplpy.core]
INFO: Auto-setting vmax to  3.771e+00 [aplpy.core]
Moment 0 for sample pixel Region 6 is 1.3153144890763757 Jy km / (beam s)
WARNING: AstropyDeprecationWarning: CoordinateHelper.ticks should not be accessed directly and is deprecated [astropy.visualization.wcsaxes.coordinate_helpers]
INFO: Auto-setting vmin to -3.355e+03 [aplpy.core]
INFO: Auto-setting vmax to  4.894e+03 [aplpy.core]
Moment 1 for sample pixel Region 6 is 178.77670374261976 km / s
WARNING: AstropyDeprecationWarning: CoordinateHelper.ticks should not be accessed directly and is deprecated [astropy.visualization.wcsaxes.coordinate_helpers]
INFO: Auto-setting vmin to -2.090e+01 [aplpy.core]
INFO: Auto-setting vmax to  3.980e+02 [aplpy.core]
Moment 2 for sample pixel Region 6 is 103.65347329551717 km / s
```

As you can see, **CubeLineMoment** is very chatty.  We will likely cut this back this verbosity a bit at some point in the future.  Note that all of the warnings are ignorable (see [Masking Used in CubeLineMoment](#masking-used-in-cubelinemoment)), due to minor things like using NaNs for blanking in the input image cube.  What you should see how are a number of new directories: 

```sh
% ls *.png
DEBUG_plot_NGC253_HCN_32_widthscale1.0_sncut999.0_widthcutscale1.0.png
```

...and there should be five new directories...
```
diagnostics
moment0
moment1
moment2
subcubes
```
The <span style="color:blue">diagnostics, moment0, moment1, moment2, and subcubes</span> directories have been created by **CubeLineMoment**.  These directories contain:
* **diagnostics:** Diagnostic plots of the spectrum extraction for each transition requested in the **yaml* input file which shows the masking used and an example gaussian fit to the respective transition.  For example, the diagnostic plot for the GMC6 position from the sample pixel regions file for the HCN 1-0 transition looks like this: ![](./images/NGC253_HCN_10_widthscale1.0_sncut3.0_widthcutscale1.0_spectraldiagnostics_GMC6.png "HCN 1-0 diagnotics example")  Note how the gaussian fit to this transition would not have been very good.  A PPV masking diagnostic plot is also produced which looks like the following for the current example: ![](./images/NGC253_brightest_diagnostic_samplepixel_GMC6.png "PPV diagnostic example")
* **moment0:** The derived zeroth moment (integrated intensity) images for each transition in FITS and png format.  The moment0 image from the current example looks like this: ![](./images/NGC253_HCN_10_moment0_widthscale1.0_sncut3.0_widthcutscale1.0.png "HCN 1-0 moment0") This directory also contains several diagnostic FITS files which show the various masking parameters used: **CentroidMap**, **FWHMMap**, **MaxMap**, **NoiseMap**, **NoiseMapBright**, and **WidthMap**.
* **moment1:** The derived first moment (centroid velocity) images for each transition in FITS and png format.  The moment1 image from the current example looks like this: ![](./images/NGC253_HCN_10_moment1_widthscale1.0_sncut3.0_widthcutscale1.0.png "HCN 1-0 moment1")
* **moment2:** The derived second moment (velocity width) images for each transition in FITS and png format.  The moment2 image from the current example looks like this: ![](./images/NGC253_HCN_10_moment2_widthscale1.0_sncut3.0_widthcutscale1.0.png "HCN 1-0 moment2")
* **subcubes:** The derived subcubes for each transition derived using the specified masking in FITS format.  In other words, these are single-transition spectral line cubes of all transitions requested in the input **yaml** file.


```python

```
