"""

Derive Moment0, Moment1, and Moment2 from a reasonably-well separated spectral line in
an image cube.  Simply calculates moments over a defined HWZI for each line in band.

To run in ipython use:

>>> %run CubeLineMoment.py inputfile.yaml

"""
from __future__ import print_function

import os
import numpy as np
from spectral_cube import SpectralCube
from astropy import units as u
from astropy import constants
import regions
import pylab as pl
import yaml
import warnings
import ast
from astropy import wcs

warnings.filterwarnings('ignore', category=wcs.FITSFixedWarning)

#from astropy import log
#log.setLevel('CRITICAL') # disable most logger messages

# suppress errors when we do np.nan > 5, etc.
np.seterr(invalid='ignore')


# debugging:
from astropy import log
import psutil
proc = psutil.Process()



def cubelinemoment_setup(cube, cuberegion, cutoutcube,
                         cutoutcuberegion, vz, target, brightest_line_frequency,
                         width_line_frequency, velocity_half_range,
                         noisemapbright_baseline, noisemap_baseline,
                         spatial_mask_limit, mask_negatives=True,
                         sample_pixel=None, **kwargs):
    """
    For a given cube file, read it and compute the moments (0,1,2) for a
    selection of spectral lines.  This code is highly configurable.

    In the parameter description, 'PPV' refers to position-position-velocity,
    and all cubes are expected to be in this space.  Velocity is
    generally interchangeable with frequency, but many operations must be
    performed in velocity space.

    Parameters
    ----------
    cube : str
        The cube file name
    cuberegion : str, optional
        A ds9 region file specifying a spatial region to extract from the cube
    cutoutcube : str
        Filename of a cube that specifies the PPV region over which the moments
        will be extracted.
    cutoutcuberegion : str, optional
        A ds9 region file specifying a spatial region to extract from the
        spatial mask cube.  Should generally be the same as cuberegion.
        NOTE TO JEFF: should this *always* be the same as cuberegion?
    vz : `astropy.units.Quantity` with km/s equivalence
        The line-of-sight velocity of the source, e.g., the redshift.
    target : str
        Name of the source.  Used when writing output files.
    brightest_line_frequency : `astropy.units.Quantity` with Hz equivalence
        The frequency of the brightest line, used to establish the cube volume
        over which to compute moments for other lines
    width_line_frequency : `astropy.units.Quantity` with Hz equivalence
        The central frequency of the line used to compute the width (moment 2)
    velocity_half_range : `astropy.units.Quantity` with km/s equivalence
        The approximate half-width zero-intensity of the lines.  This parameter
        is used to crop out regions of the cubes around line centers.  It
        should be larger than the expected FWHM line width.  It should encompass
        the *full width* of the line over the *whole* source, i.e., if your
        Galaxy has a rotation curve from -100 to +100 km/s and a typical LOS
        linewidth of 20 km/s, it should go from -120 to +120 (or may be -140 to
        +140 to be conservative)
    noisemapbright_baseline : list of lists
        A list of pairs of indices over which the noise can be computed from
        the 'bright' cube
        NOTE TO JEFF: It would probably be better to specify this in GHz or
        km/s.  That will require a slight change in the code, but will make
        it more robust to changes in, e.g., linewidth or other parameters
        that can affect the cube shape.
    noisemap_baseline : list of lists
        A list of pairs of indices over which the noise can be computed from
        the main cube
    spatial_mask_limit : float
        Factor in n-sigma above which to apply threshold.  Any spatial pixels
        whose *peak intensity* is below this limit will be flagged out.
    mask_negatives : float or bool
        Mask out negatives below N-sigma negative.
    sample_pixel : tuple
        A set of (x,y) coordinates to sample from the cutout cube to create
        diagnostic images.  Note that these must be *in the frame of the cutout
        cube*, which means you may need to do a little math to make sure
        they're correct.  If left as `None`, no diagnostic images will be made.


    Returns
    -------
    A variety of cubes and maps
    """

    # Read the FITS cube
    # And change the units back to Hz
    cube = SpectralCube.read(cube).with_spectral_unit(u.Hz)

    # cut out a region that only includes the Galaxy (so we don't have to worry
    # about masking later)
    if cuberegion is not None:
        cube = cube.subcube_from_regions(regions.read_ds9(cuberegion))

    # --------------------------
    # Define a spatial mask that guides later calculations by defining where
    # dense gas is and is not.
    # For the NGC253 Band 6 data use the C18O 2-1 line in spw1 for the dense
    # gas mask for all Band 6 lines.
    #    cutoutcube = SpectralCube.read('NGC253-H213COJ32K1-Feather-line-All.fits').with_spectral_unit(u.Hz).subcube_from_regions(regions.read_ds9('ngc253boxband6tight.reg'))
    cutoutcube = (SpectralCube.read(cutoutcube)
                  .with_spectral_unit(u.Hz)
                  )
    if cutoutcuberegion is not None:
        cutoutcube = cutoutcube.subcube_from_regions(regions.read_ds9(cutoutcuberegion))
    noisecube = cutoutcube
    # For the NGC4945 Band 6 data use the C18O 2-1 line in spw1 for the dense
    # gas mask for all Band 6 lines.
    #cutoutcube = SpectralCube.read('NGC4945-H213COJ32K1-Feather-line.fits').with_spectral_unit(u.Hz).subcube_from_regions(regions.read_ds9('ngc4945boxband6.reg'))

    if mask_negatives is not False:
        std = cube.std()
        posmask = cutoutcube > (std * mask_negatives)
        cutoutcube = cutoutcube.with_mask(posmask)


    # redshift velocity
    #    vz = 258.8*u.km/u.s # For NGC253
    vz = u.Quantity(vz, u.km/u.s) # For NGC253
    #vz = 538.2*u.km/u.s # For NGC4945

    # Lines to be analyzed (including brightest_line)
    #    target = 'NGC253'
    #target = 'NGC4945'

    #    brightest_line_frequency = 219.560358*u.GHz # C18O 2-1
    brightest_line_frequency = u.Quantity(brightest_line_frequency, u.GHz) # C18O 2-1
    #    width_line = 218.222192*u.GHz # H2CO 3(03)-2(02)
    # NOT USED width_line_frequency = u.Quantity(width_line_frequency, u.GHz) # H2CO 3(03)-2(02)

    # Assume you have a constant expected width (HWZI) for the brightest line
    # Note: This HWZI should be larger than those assumed in the line extraction loop below...
    #    width = 80*u.km/u.s
    velocity_half_range = u.Quantity(velocity_half_range, u.km/u.s)

    # Create a copy of the cutoutcube with velocity units
    cutoutVcube = cutoutcube.with_spectral_unit(u.km/u.s,
                                                   rest_value=brightest_line_frequency,
                                                   velocity_convention='optical')

    # Use the brightest line to identify the appropriate peak velocities, but ONLY
    # from a slab including +/- width:
    brightest_cube = cutoutVcube.spectral_slab(vz-velocity_half_range,
                                               vz+velocity_half_range)

    # compute various moments & statistics along the spcetral dimension
    peak_velocity = brightest_cube.spectral_axis[brightest_cube.argmax(axis=0)]
    max_map = peak_amplitude = brightest_cube.max(axis=0)
    width_map = brightest_cube.linewidth_sigma() # or vcube.moment2(axis=0)**0.5
    fwhm_map = brightest_cube.linewidth_fwhm() # FOR TESTING
    sqrtmom2_map = brightest_cube.moment2(axis=0)**0.5 # FOR TESTING
    centroid_map = brightest_cube.moment1(axis=0)

    # NOTE: the updating header stuff will be completely redundant after
    # https://github.com/radio-astro-tools/spectral-cube/pull/383 is merged
    if not os.path.exists('moment0'):
        os.mkdir('moment0')
    hdu = width_map.hdu
    hdu.header['OBJECT'] = cube.header['OBJECT']
    hdu.writeto("moment0/{0}_WidthMap.fits".format(target),overwrite=True)
    hdu = centroid_map.hdu
    hdu.header['OBJECT'] = cube.header['OBJECT']
    hdu.writeto("moment0/{0}_CentroidMap.fits".format(target),overwrite=True)
    hdu = peak_amplitude.hdu
    hdu.header['OBJECT'] = cube.header['OBJECT']
    hdu.writeto("moment0/{0}_MaxMap.fits".format(target),overwrite=True)
    hdu = fwhm_map.hdu
    hdu.header['OBJECT'] = cube.header['OBJECT']
    hdu.writeto("moment0/{0}_FWHMMap.fits".format(target),overwrite=True)
    hdu = sqrtmom2_map.hdu
    hdu.header['OBJECT'] = cube.header['OBJECT']
    hdu.writeto("moment0/{0}_SQRTMOM2Map.fits".format(target),overwrite=True)


    # From NGC253 H213COJ32K1 spectral baseline
    inds = np.arange(noisecube.shape[0])
    mask = np.zeros_like(inds, dtype='bool')
    for low,high in noisemapbright_baseline:
        mask[low:high] = True

    # need to use an unmasked cube
    noisemapbright = noisecube.with_mask(mask[:,None,None]).std(axis=0)
    print("noisemapbright peak = {0}".format(np.nanmax(noisemapbright)))

    # Make a plot of the noise map...
    #pl.figure(2).clf()
    #pl.imshow(noisemapbright.value)
    #pl.colorbar()
    hdu = noisemapbright.hdu
    hdu.header.update(cutoutcube.beam.to_header_keywords())
    hdu.header['OBJECT'] = cutoutcube.header['OBJECT']
    hdu.writeto("moment0/{0}_NoiseMapBright.fits".format(target),overwrite=True)
    #
    # Use 3*noisemap for spatial masking
    if spatial_mask_limit is None:
        spatial_mask = np.ones(noisemapbright.shape, dtype='bool')
    else:
        spatial_mask = np.fabs(peak_amplitude) > spatial_mask_limit*noisemapbright
    #hdu = spatial_mask.hdu
    #hdu.header.update(cutoutcube.beam.to_header_keywords())
    #hdu.header['OBJECT'] = cutoutcube.header['OBJECT']
    #hdu.writeto("moment0/{0}_ppvmask.fits".format(target),overwrite=True)
    # --------------------------

    # Now process spw of interest...
    #
    # Now define noise map for spw being analyzed...
    # From NGC253 H2COJ32K02 spectral baseline
    #noisemap = cube[360:370,:,:].std(axis=0)
    # ADAM ADDED: Derive noisemap over non-contiguous baseline
    # JGM: Had to go back to defining noisemap_baseline in function as param input of list does not seem to work
    #noisemap_baseline = [(9, 14), (40, 42), (72, 74), (114, 122), (138, 143), (245, 254), (342, 364)]
    inds = np.arange(cube.shape[0])
    mask = np.zeros_like(inds, dtype='bool')
    for low,high in noisemap_baseline:
        mask[low:high] = True
    noisemap = cube.with_mask(mask[:,None,None]).std(axis=0)
    hdu = noisemap.hdu
    hdu.header.update(cube.beam.to_header_keywords())
    hdu.header['OBJECT'] = cube.header['OBJECT']
    hdu.writeto("moment0/{0}_NoiseMap.fits".format(target),overwrite=True)

    if sample_pixel is not None:
        # Create a plot showing all the analysis steps applied to the sample
        # pixel
        fig = pl.figure(11)
        fig.clf()
        #ax = fig.gca() # subplot?
        #raw_spec = cube[:, sample_pixel[0], sample_pixel[1]]
        #ax.plot(raw_spec.spectral_axis, raw_spec.value, drawstyle='steps-mid',
        #        color='k', label='Raw')
        ax = fig.add_subplot(3,1,1)
        ppvmaskplot = cutoutcube[:, sample_pixel[0], sample_pixel[1]]
        ax.plot(ppvmaskplot.spectral_axis, ppvmaskplot.value,
                drawstyle='steps-mid', color='k', label='PPV Mask')
        ax.set_title('PPV Mask')

        ax2 = fig.add_subplot(3,1,2)
        noisespec = noisecube[:, sample_pixel[0], sample_pixel[1]]
        ax2.plot(noisespec.spectral_axis, noisespec.value,
                 drawstyle='steps-mid', color='b', label='Noise')
        ax2.set_title('Noise')

        ax3 = fig.add_subplot(3,1,3)
        brightestspec = brightest_cube[:, sample_pixel[0], sample_pixel[1]]
        ax3.plot(brightestspec.spectral_axis, brightestspec.value,
                 drawstyle='steps-mid', color='r', label='Brightest')
        ax3.set_title('Brightest')

        ax2.plot(brightest_cube.with_spectral_unit(cutoutcube.spectral_axis.unit).spectral_axis.value,
                 brightestspec.value,
                 drawstyle='steps-mid', color='r', label='Brightest',
                 zorder=-1, linewidth=2)

        if not os.path.exists('diagnostics'):
            os.mkdir('diagnostics')

        fig.savefig('diagnostics/{0}_brightest_diagnostic.png'.format(target))

    return (cube, cutoutcube, spatial_mask, noisemap, noisemapbright,
            centroid_map, width_map, max_map, peak_velocity)



def cubelinemoment_multiline(cube, peak_velocity, centroid_map, max_map,
                             noisemap, signal_mask_limit,
                             my_line_list, my_line_widths, my_line_names,
                             target, spatial_mask, width_map,
                             width_map_scaling=1.0, width_cut_scaling=1.0,
                             fit=False, apply_width_mask=True,
                             sample_pixel=None,
                             **kwargs):
    """
    Given the appropriate setup, extract moment maps for each of the specified
    lines

    Parameters
    ----------
    peak_velocity : `astropy.units.Quantity` with km/s equivalence
    centroid_map : `astropy.units.Quantity` with km/s equivalence
    max_map : `astropy.units.Quantity` with brightness or flux unit
    noisemap : `astropy.units.Quantity` with brightness or flux unit
    my_line_list : `astropy.units.Quantity` with Hz equivalence
        An array of line centers to compute the moments of
    my_line_widths : `astropy.units.Quantity` with km/s equivalence
        An array of line widths matched to ``my_line_list``.
    my_line_names : list of strings
        A list of names matched to ``my_line_list`` and ``my_line_widths``.
        Used to specify the output filename.
    signal_mask_limit : float
        Factor in n-sigma above which to apply threshold to data.  Unlike
        ``spatial_mask_limit``, this threshold is applied on a per-voxel basis.
        If this is set to ``None``, no signal masking will be applied.
    width_map_scaling : float
        A factor by which to multiply the ``width_map`` when making the
        position-velocity mask cube.
    width_cut_scaling : float
        The factor by which the cube cutout is expanded, so if this is != 1,
        the extracted subcube will be larger.
    apply_width_mask : bool
        Should width masking be applied at all?  Turning this off can save some
        computational time.

    Returns
    -------
    None.  Outputs are saved to files in the momentX/ subdirectory,
    where X is in {0,1,2}
    """

    # parameter checking
    if len(my_line_names) != len(my_line_list) or len(my_line_names) != len(my_line_widths):
        raise ValueError("Line lists (central frequency, names, and widths) "
                         "have different lengths")


    # Now loop over EACH line, extracting moments etc. from the appropriate region:
    # we'll also apply a transition-dependent width (my_line_widths) here because
    # these fainter lines do not have peaks as far out as the bright line.

    for line_name,line_freq,line_width in zip(my_line_names,my_line_list,my_line_widths):

        log.info("Line: {0}, {1}, {2}".format(line_name, line_freq, line_width))

        line_freq = u.Quantity(line_freq,u.GHz)
        line_width = u.Quantity(line_width,u.km/u.s) * width_cut_scaling
        vcube = cube.with_spectral_unit(u.km/u.s, rest_value=line_freq,
                                        velocity_convention='optical')

        subcube = vcube.spectral_slab(peak_velocity.min()-line_width,
                                      peak_velocity.max()+line_width)

        if apply_width_mask:
            # ADAM'S ADDITIONS AGAIN
            # use the spectral_axis to make a 'mask cube' with the moment1/moment2
            # values computed for the selected mask line (H2CO 303?)
            # We create a Gaussian along each line-of-sight, then we'll crop based on a
            # threshold
            # The [:,:,None] and [None,None,:] allow arrays of shape [x,y,0] and
            # [0,0,z] to be "broadcast" together
            assert centroid_map.unit.is_equivalent(u.km/u.s)
            gauss_mask_cube = np.exp(-(np.array(centroid_map)[None,:,:] -
                                       np.array(subcube.spectral_axis)[:,None,None])**2 /
                                     (2*np.array(width_map*width_map_scaling)[None,:,:]**2))
            peak_sn = max_map / noisemap

            print("Peak S/N: {0}".format(np.nanmax(peak_sn)))

            # threshold at the fraction of the Gaussian corresponding to our peak s/n.
            # i.e., if the S/N=6, then the threshold will be 6-sigma
            # (this can be modified as you see fit)
            threshold = 1 / peak_sn

            print("Highest Threshold: {0}".format(np.nanmax(threshold)))
            #print("Lowest Threshold: {0}".format((threshold[threshold>0].min())))
            if sample_pixel:
                print("SP Threshold: {0}".format(threshold[sample_pixel]))
                print("SP S, N, S/N: {0},{1},{2}"
                      .format(max_map[sample_pixel],
                              noisemap[sample_pixel],
                              peak_sn[sample_pixel],
                             ))


            # this will compare the gaussian cube to the threshold on a (spatial)
            # pixel-by-pixel basis
            width_mask_cube = gauss_mask_cube > threshold
            print("Number of values above threshold: {0}".format(width_mask_cube.sum()))
            print("Max value in the mask cube: {0}".format(np.nanmax(gauss_mask_cube)))
            print("shapes: mask cube={0}  threshold: {1}".format(gauss_mask_cube.shape, threshold.shape))

            msubcube = subcube.with_mask(width_mask_cube)
        else:
            msubcube = subcube


        # Mask on a pixel-by-pixel basis with an N-sigma cut
        if signal_mask_limit is not None:
            signal_mask = subcube > signal_mask_limit*noisemap
            msubcube = msubcube.with_mask(signal_mask)


        # this part makes a cube of velocities
        temp = subcube.spectral_axis
        velocities = np.tile(temp[:,None,None], subcube.shape[1:])

        # now we use the velocities from the brightest line to create a mask region
        # in the same velocity range but with different rest frequencies (different
        # lines)
        velocity_range_mask = np.abs(peak_velocity - velocities) < line_width
        # the mask is a cube, the spatial mask is a 2d array, but in this case
        # numpy knows how to combine them properly
        # (signal_mask is a different type, so it can't be combined with the others
        # yet - I'll add a feature request for that)
        msubcube = msubcube.with_mask(velocity_range_mask & spatial_mask)

        # DEBUG: show the values from all the masks
        pl.figure(10).clf()
        pl.subplot(2,2,1).imshow(velocity_range_mask.max(axis=0), origin='lower', interpolation='nearest')
        pl.subplot(2,2,1).set_title("velocity range mask")
        pl.subplot(2,2,2).imshow(spatial_mask, origin='lower', interpolation='nearest')
        pl.subplot(2,2,2).set_title("spatial mask")
        if signal_mask_limit is not None:
            pl.subplot(2,2,3).imshow(signal_mask.include().max(axis=0), origin='lower', interpolation='nearest')
        pl.subplot(2,2,3).set_title("signal mask")
        if apply_width_mask:
            pl.subplot(2,2,4).imshow(width_mask_cube.max(axis=0), origin='lower', interpolation='nearest')
        pl.subplot(2,2,4).set_title("width mask")
        pl.savefig("DEBUG_plot_{0}_{1}_widthscale{2:0.1f}_sncut{3:0.1f}_widthcutscale{4:0.1f}.png"
                   .format(target, line_name, width_map_scaling,
                           signal_mask_limit or 999, width_cut_scaling))

        if sample_pixel is not None:
            # Create a plot showing all the analysis steps applied to the sample
            # pixel
            fig = pl.figure(11)
            fig.clf()
            #ax = fig.gca() # subplot?
            #raw_spec = cube[:, sample_pixel[0], sample_pixel[1]]
            #ax.plot(raw_spec.spectral_axis, raw_spec.value, drawstyle='steps-mid',
            #        color='k', label='Raw')
            ax1 = fig.add_subplot(2,1,1)
            subcubesp = subcube[:, sample_pixel[0], sample_pixel[1]]
            ax1.plot(subcubesp.spectral_axis, subcubesp.value,
                     drawstyle='steps-mid', color='k', label='subcube')
            ax1.set_title('subcube')

            ax = fig.add_subplot(2,1,2)
            mask_ = msubcube.mask.include()[:, sample_pixel[0], sample_pixel[1]]
            maskedsubcubesp = msubcube[:, sample_pixel[0], sample_pixel[1]]
            assert np.all(np.isfinite(maskedsubcubesp[mask_]))
            assert np.all(~np.isfinite(maskedsubcubesp[~mask_]))
            nansp = maskedsubcubesp.filled_data[:]
            zerosp = np.nan_to_num(nansp)
            ax.plot(subcubesp.spectral_axis, subcubesp.value,
                    drawstyle='steps-mid', linestyle=":", color='k',
                    zorder=-30,
                    label='subcube')
            ax.plot(maskedsubcubesp.spectral_axis, nansp.value,
                    linewidth=7, zorder=-50,
                    drawstyle='steps-mid', color='k', alpha=0.3,
                    label='Masked subcube')
            ax.set_title('masked subcube')

            ax.plot(velocities[:, sample_pixel[0], sample_pixel[1]],
                    subcubesp.value*velocity_range_mask[:, sample_pixel[0], sample_pixel[1]],
                    color='orange',
                    linewidth=3,
                    zorder=-15,
                    alpha=0.5,
                    label='VelocityRangeMask',
                    drawstyle='steps-mid',
                   )


            if apply_width_mask and 'width_mask_cube' in locals():
                ax.plot(maskedsubcubesp.spectral_axis,
                        subcubesp.value*width_mask_cube[:, sample_pixel[0], sample_pixel[1]],
                        drawstyle='steps-mid', color='b', label='Width Mask',
                        alpha=0.5, zorder=-10, linewidth=3)
                ax.plot(maskedsubcubesp.spectral_axis,
                        gauss_mask_cube[:, sample_pixel[0], sample_pixel[1]] * subcubesp.value.max(),
                        color='r', zorder=-20, linewidth=1,
                        label='Gaussian',
                       )
            if 'signal_mask' in locals():
                ax.plot(maskedsubcubesp.spectral_axis,
                        subcubesp.value*signal_mask[:, sample_pixel[0], sample_pixel[1]].include(),
                        drawstyle='steps-mid', color='g', label='Signal Mask',
                        alpha=0.5, zorder=-10, linewidth=3)

            pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))

            if not os.path.exists('diagnostics'):
                os.mkdir('diagnostics')

            fig.savefig("diagnostics/{0}_{1}_widthscale{2:0.1f}_sncut{3:0.1f}_widthcutscale{4:0.1f}_spectraldiagnostics.png"
                        .format(target, line_name, width_map_scaling,
                                signal_mask_limit or 999, width_cut_scaling),
                        bbox_inches='tight')

        # Now write output.  Note that moment0, moment1, and moment2 directories
        # must already exist...

        labels = {0: 'Integrated Intensity [{0}]',
                  1: '$V_{{LSR}}$ [{0}]',
                  #2: '$\sigma_v$ [{0}]',
                  2: '$FWHM$ [{0}]',
                 }

        moments = {}

        pl.close('all')

        for moment in (0,1,2):
            if not os.path.exists('moment{0}'.format(moment)):
                os.mkdir('moment{0}'.format(moment))
            if moment == 2:
                mom = msubcube.linewidth_fwhm()
            else:
                mom = msubcube.moment(order=moment, axis=0)
            hdu = mom.hdu
            hdu.header.update(cube.beam.to_header_keywords())
            hdu.header['OBJECT'] = cube.header['OBJECT']
            hdu.writeto("moment{0}/{1}_{2}_moment{0}_widthscale{3:0.1f}_sncut{4:0.1f}_widthcutscale{5:0.1f}.fits"
                        .format(moment, target, line_name, width_map_scaling,
                                signal_mask_limit or 999, width_cut_scaling), overwrite=True)
            pl.figure(1).clf()
            mom.quicklook() #filename='moment{0}/{1}_{2}_moment{0}.png'.format(moment,target,line_name))
            mom.FITSFigure.colorbar.show(axis_label_text=labels[moment].format(mom.unit.to_string('latex_inline')))
            mom.FITSFigure.save(filename='moment{0}/{1}_{2}_moment{0}_widthscale{3:0.1f}_sncut{4:0.1f}_widthcutscale{5:0.1f}.png'
                                .format(moment, target, line_name,
                                        width_map_scaling, signal_mask_limit or 999, width_cut_scaling))
            mom.FITSFigure.close()
            moments[moment] = mom

            if sample_pixel is not None:
                print("Moment {0} for sample pixel is {1}"
                      .format(moment, mom[sample_pixel]))

        if not os.path.exists('subcubes'):
            os.mkdir('subcubes')

        subcube_outname = ('subcubes/{0}_{1}_widthscale{4:0.1f}_widthcutscale{2:0.1f}_sncut{3:0.1f}_subcube.fits'
                           .format(target, line_name, width_cut_scaling,
                                   signal_mask_limit or 999, width_map_scaling))
        msubcube.write(subcube_outname, overwrite=True)

        # finally, optionally, do some pyspeckit fitting
        if fit:
            import pyspeckit
            msubcube_allvalid = msubcube._new_cube_with()
            msubcube_allvalid._mask = None
            pcube = pyspeckit.Cube(cube=msubcube)
            max_map_sub = msubcube.max(axis=0).value
            pcube.mapplot.plane = max_map_sub
            guesses = np.array([max_map_sub, moments[1].value,
                                moments[2].value / (8*np.log(2))**0.5])
            maskmap = (np.all(guesses > 0, axis=0) &
                       (msubcube.mask.include().sum(axis=0) > 3))
            print("Fitting {0} spectra with pyspeckit".format(maskmap.sum()))
            pcube.fiteach(guesses=guesses, start_from_point='center',
                          errmap=noisemap.value, signal_cut=0, maskmap=maskmap,
                          limited=[(True,True),(True,True),(True,True)],
                          limits=[(0,max_map_sub.max()*2),
                                  (moments[1].value.min()-50, moments[1].value.max()+50),
                                  (0, guesses[2,:,:].max()*2)],
                         )
            pcube.write_fit('pyspeckit_fits/{0}_{1}_fitcube.fits'.format(target,
                                                                         line_name),
                            overwrite=True)
        log.debug("Open files: {0}".format(len(proc.open_files())))

    return locals()


def pyspeckit_fit_cube(cube, max_map, centroid_map, width_map, noisemap,
                       lines, vz):
    """
    This is experimental and doesn't really work: the idea here is to fit all
    lines in the cube simultaneously.
    """
    import pyspeckit

    vz = u.Quantity(vz, u.km/u.s)

    fcube = cube.with_spectral_unit(u.GHz)

    def inrange(x):
        return (x < fcube.spectral_extrema[1] and
                x > fcube.spectral_extrema[0])

    lines_in_cube = {linename: linedata
                     for linename, linedata in lines.items()
                     if inrange(linedata['frequency']*(1-vz/constants.c))}

    frequencies = sorted(linedata['frequency'] for linedata in lines_in_cube.values())

    line_guesses = [[max_map.value,
                     ((1-centroid_map/constants.c)*frq).to(u.GHz).value,
                     ((width_map/constants.c)*frq).to(u.GHz).value]
                    for frq in frequencies]
    line_guesses = np.array([x for y in line_guesses for x in y])

    guesses = np.array([max_map.value, centroid_map.value, width_map.value])
    #vcube = cube.with_spectral_unit(u.km/u.s, velocity_convention='optical')
    pcube = pyspeckit.Cube(cube=fcube)
    pcube.mapplot.plane = max_map.value
    pcube.fiteach(guesses=guesses, start_from_point=(150,150),
                  errmap=noisemap.value)


def main():
    """
    To avoid ridiculous namespace clashes
    http://stackoverflow.com/questions/4775579/main-and-scoping-in-python
    """

    import argparse

    parser = argparse.ArgumentParser(description='Derive moment maps for a'
                                     ' cube given a complex suite of'
                                     ' parameters')
    parser.add_argument('param_file', metavar='pars', type=str,
                        help='The name of the YAML parameter file')

    args = parser.parse_args()

    infile = args.param_file

    # Read input file which sets all parameters for processing
    # Example call:
    # ipython:
    # %run CubeLineMoment.py yaml_scripts/NGC253-H2COJ32K02-CubeLineMomentInput.yaml
    # cmdline:
    # python CubeLineMoment.py yaml_scripts/NGC253-H2COJ32K02-CubeLineMomentInput.yaml

    with open(infile) as fh:
        params = yaml.load(fh)



    for par in params:
        if params[par] == 'None':
            params[par] = None

    if params['signal_mask_limit'] == 'None':
        params['signal_mask_limit'] = None
    elif hasattr(params['signal_mask_limit'], 'split'):
        params['signal_mask_limit'] = list(map(float, params['signal_mask_limit'].split(", ")))
    if params['spatial_mask_limit'] == 'None':
        params['spatial_mask_limit'] = None
    elif hasattr(params['spatial_mask_limit'], 'split'):
        params['spatial_mask_limit'] = list(map(float, params['spatial_mask_limit'].split(", ")))
    if 'width_map_scaling' in params and hasattr(params['width_map_scaling'], 'split'):
        params['width_map_scaling'] = list(map(float, params['width_map_scaling'].split(", ")))
    if 'width_cut_scaling' in params and hasattr(params['width_cut_scaling'], 'split'):
        params['width_cut_scaling'] = list(map(float, params['width_cut_scaling'].split(", ")))
    params['my_line_list'] = u.Quantity(list(map(float, params['my_line_list'].split(", "))), u.GHz)
    params['my_line_widths'] = u.Quantity(list(map(float, params['my_line_widths'].split(", "))), u.km/u.s)
    params['my_line_names'] = params['my_line_names'].split(", ")
    if 'sample_pixel' in params:
        params['sample_pixel'] = ast.literal_eval(params['sample_pixel'])

    print(params)

    # Read parameters from dictionary

    (cube, spatialmaskcube, spatial_mask, noisemap, noisemapbright,
     centroid_map, width_map, max_map, peak_velocity) = cubelinemoment_setup(**params)

    params.pop('cube')

    cubelinemoment_multiline(cube=cube, spatial_mask=spatial_mask,
                             peak_velocity=peak_velocity,
                             centroid_map=centroid_map, max_map=max_map,
                             noisemap=noisemap, width_map=width_map, fit=False,
                             **params)

    # params.pop('signal_mask_limit')
    # cubelinemoment_multiline(cube=cube, spatial_mask=spatial_mask,
    #                          peak_velocity=peak_velocity,
    #                          centroid_map=centroid_map, max_map=max_map,
    #                          noisemap=noisemap, width_map=width_map,
    #                          width_map_scaling=2.0, fit=False,
    #                          signal_mask_limit=2.0,
    #                          **params)

    # cubelinemoment_multiline(cube=cube, spatial_mask=spatial_mask,
    #                          peak_velocity=peak_velocity,
    #                          centroid_map=centroid_map, max_map=max_map,
    #                          noisemap=noisemap, width_map=width_map,
    #                          width_map_scaling=2.0, fit=False,
    #                          width_cut_scaling=1.5,
    #                          signal_mask_limit=2.0,
    #                          **params)

    # cubelinemoment_multiline(cube=cube, spatial_mask=spatial_mask,
    #                          peak_velocity=peak_velocity,
    #                          centroid_map=centroid_map, max_map=max_map,
    #                          noisemap=noisemap, width_map=width_map,
    #                          width_map_scaling=1.0, fit=False,
    #                          width_cut_scaling=1.0,
    #                          signal_mask_limit=2.0,
    #                          **params)

    # Clean up open figures
    pl.close('all')

    # useful reformatting of the lines to pass to the pyspeckit fitter if we
    # ever choose to use it
    lines = dict(zip(params['my_line_names'],
                     [{'frequency':frq,
                      'width':wid}
                      for frq,wid in zip(params['my_line_list'],
                                         params['my_line_widths'])]
                    ))

    return locals()


if __name__ == "__main__":
    new_locals = main()
    locals().update(new_locals)

