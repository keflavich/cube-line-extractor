"""

Derive Moment0, Moment1, and Moment2 from a reasonably-well separated spectral line in
an image cube.  Simply calculates moments over a defined HWZI for each line in band.

Band 6 version (also have a Band 7 specific version)

To run in ipython use:

run ~/Python/CubeLineMoment.py

"""
from __future__ import print_function

import numpy as np
from spectral_cube import SpectralCube
from astropy import wcs
from astropy import units as u
import pyregion
import pylab as pl
import aplpy
import radio_beam

from astropy import log
log.setLevel('CRITICAL') # disable most logger messages

import yaml

# Read input file which sets all parameters for processing
infile = raw_input('Enter input file name: ')
with open(infile) as fh:
    params = yaml.load(fh)

print(params)

def cubelinemoment(cube, cuberegion, spatialmaskcube, spatialmaskcuberegion,
                   vz, target, brightest_line, width_line, width,
                   noisemapbright_baseline, noisemap_baseline, my_line_list,
                   my_line_widths, my_line_names, signal_mask_limit,
                   spatial_mask_limit):

    # Read the FITS cube
    # And change the units back to Hz
    # and cut out a region that only includes the Galaxy (so we don't have to worry
    # about masking later)
    #
    # The following for NGC253-H2COJ32K02...
    #    cube = SpectralCube.read('NGC253-H2COJ32K02-Feather-line-All.fits').with_spectral_unit(u.Hz).subcube_from_ds9region(pyregion.open('ngc253boxband6tight.reg'))
    cube = SpectralCube.read(cube).with_spectral_unit(u.Hz).subcube_from_ds9region(pyregion.open(cuberegion))
    # The following for NGC4945-H2COJ32K02...
    #cube = SpectralCube.read('NGC4945-H2COJ32K02-Feather-line.fits').with_spectral_unit(u.Hz).subcube_from_ds9region(pyregion.open('ngc4945boxband6.reg'))
    # The following for NGC253-H2COJ54K23...
    #cube = SpectralCube.read('NGC253-H2COJ54K23-Feather-line.fits').with_spectral_unit(u.Hz).subcube_from_ds9region(pyregion.open('ngc253box.reg'))
    # The following for NGC253-H2COJ54K1...
    #cube = SpectralCube.read('NGC253-H2COJ54K1-Feather-line.fits').with_spectral_unit(u.Hz).subcube_from_ds9region(pyregion.open('ngc253box.reg'))

    # --------------------------
    # Define a spatial mask that guides later calculations by defining where
    # dense gas is and is not.
    # For the NGC253 Band 6 data use the C18O 2-1 line in spw1 for the dense
    # gas mask for all Band 6 lines.
    #    spatialmaskcube = SpectralCube.read('NGC253-H213COJ32K1-Feather-line-All.fits').with_spectral_unit(u.Hz).subcube_from_ds9region(pyregion.open('ngc253boxband6tight.reg'))
    spatialmaskcube = SpectralCube.read(spatialmaskcube).with_spectral_unit(u.Hz).subcube_from_ds9region(pyregion.open(spatialmaskcuberegion))
    # For the NGC4945 Band 6 data use the C18O 2-1 line in spw1 for the dense
    # gas mask for all Band 6 lines.
    #spatialmaskcube = SpectralCube.read('NGC4945-H213COJ32K1-Feather-line.fits').with_spectral_unit(u.Hz).subcube_from_ds9region(pyregion.open('ngc4945boxband6.reg'))

    # redshift velocity
    #    vz = 258.8*u.km/u.s # For NGC253
    vz = vz*u.km/u.s # For NGC253
    #vz = 538.2*u.km/u.s # For NGC4945

    # Lines to be analyzed (including brightest_line)
    #    target = 'NGC253'
    #target = 'NGC4945'

    #    brightest_line = 219.560358*u.GHz # C18O 2-1
    brightest_line = brightest_line*u.GHz # C18O 2-1
    #    width_line = 218.222192*u.GHz # H2CO 3(03)-2(02)
    width_line = width_line*u.GHz # H2CO 3(03)-2(02)

    # Assume you have a constant expected width (HWZI) for the brightest line
    # Note: This HWZI should be larger than those assumed in the line extraction loop below...
    #    width = 80*u.km/u.s
    width = width*u.km/u.s

    # ADAM'S ADDITIONS HERE
    # Use the H2CO 303_202 line (H2COJ32K02) as a mask for line widths...
    vcube = cube.with_spectral_unit(u.km/u.s, rest_value=width_line,
                                velocity_convention='optical')
    width_map = vcube.linewidth_sigma() # or vcube.moment2(axis=0)**0.5
    centroid_map = vcube.moment1(axis=0)
    max_map = cube.max(axis=0)
    #max_width = width_map.max() # should be ~150 km/s?
    #max_fwhm_width = max_width * (8*np.log(2))**0.5 # convert from sigma to FWHM


    # Use the brightest line to identify the appropriate peak velocities, but ONLY
    # from a slab including +/- width:
    brightest_cube = spatialmaskcube.with_spectral_unit(u.km/u.s,
                                                    rest_value=brightest_line,
                                                    velocity_convention='optical').spectral_slab(vz-width,
                                                                                                 vz+width)

    peak_velocity = brightest_cube.spectral_axis[brightest_cube.argmax(axis=0)]
    #pl.figure(2).clf()
    #pl.imshow(peak_velocity.value)
    #pl.colorbar()

    # make a spatial mask excluding pixels with no signal
    # (you can do better than this - this is the trivial, first try algorithm)
    peak_amplitude = brightest_cube.max(axis=0)
    # found this range from inspection of a spectrum:
    # s = cube.max(axis=(1,2))
    # s.quicklook()
    #noisemap = cube.spectral_slab(362.603*u.GHz, 363.283*u.GHz).std(axis=0)
    # Channel selection matches that used for continuum subtraction
    #
    # From NGC253 H213COJ32K1 spectral baseline
    #    noisemapbright = spatialmaskcube[150:180,:,:].std(axis=0)
    # JGM: Had to go back to defining noisemapbright_baseline in function as param input of list does not seem to work
    noisemapbright_baseline = [(150,180)]
    inds = np.arange(cube.shape[0])
    mask = np.zeros_like(inds, dtype='bool')
    for low,high in noisemapbright_baseline:
        mask[low:high] = True
    noisemapbright = cube.with_mask(mask[:,None,None]).std(axis=0)
    # From NGC4945 H213COJ32K1 spectral baseline
    #noisemapbright = spatialmaskcube[165:185,:,:].std(axis=0)
    print("noisemapbright peak = {0}".format(np.nanmax(noisemapbright)))

    # Make a plot of the noise map...
    pl.figure(2).clf()
    pl.imshow(noisemapbright.value)
    pl.colorbar()
    #
    # Use 3*noisemap for spatial masking
    spatial_mask = peak_amplitude > spatial_mask_limit*noisemapbright
    # --------------------------

    # Now process spw of interest...
    #
    # Now define noise map for spw being analyzed...
    # From NGC253 H2COJ32K02 spectral baseline
    #noisemap = cube[360:370,:,:].std(axis=0)
    # ADAM ADDED: Derive noisemap over non-contiguous baseline
    # JGM: Had to go back to defining noisemap_baseline in function as param input of list does not seem to work
    noisemap_baseline = [(9, 14), (40, 42), (72, 74), (114, 122), (138, 143), (245, 254), (342, 364)]
    inds = np.arange(cube.shape[0])
    mask = np.zeros_like(inds, dtype='bool')
    for low,high in noisemap_baseline:
        mask[low:high] = True
    noisemap = cube.with_mask(mask[:,None,None]).std(axis=0)
    #
    # From NGC4945 H2COJ32K02 spectral baseline
    #noisemap = cube[1:50,:,:].std(axis=0)
    # For NGC253 H2COJ54K23...
    #noisemap = cube[35:45,:,:].std(axis=0)
    # For NGC253 H2COJ54K1...
    #noisemap = cube[50:90,:,:].std(axis=0)
    #
    # Note that "line widths" in the following are HWZI values...
    #
    # NGC253 H2CO J=3-2 K=0,2
    #
    #    my_line_list = [217.289800, 217.299162, 217.467150, 217.517110, 217.802057, 217.88639, 217.943821, 218.15897, 218.222192, 218.324711, 218.440050, 218.475632, 218.760071, 218.85439, 218.9033555, 218.981019] * u.GHz
    #    my_line_widths = [50.0, 50.0, 60.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 50.0, 40.0, 40.0] * u.km/u.s
    #    my_line_names = ['13CNF122','CH3OH67','13CNF132','CH3OCHO88','CH3OCHO4847','CH3OH2020','CH3OCHO4546','CH3OCHO??','H2COJ32K0','HC3N2423v0','CH3OH43','H2COJ32K221','H2COJ32K210','HC3N2423v6','OCS1817','HNCO109']
    # These are:
    # 13CN N=2-1,J=5/2-3/2,F1=2-2,F=3-2 at 217289.800 MHz              (Blend with CH3OH 6(15)-7(26))
    # CH3OH 6(15)-7(26)            217299.162 MHz (+1268.05 km/s)      (Blend with 13CN F1=2-2)
    # 13CN N=2-1,J=5/2-3/2,F1=3-2,F=3-2  at 217467.150 MHz             (Set width larger to encompass multiplet)
    # CH3OCHO 8(5,4)-8(3,5)A (Methyl Formate) at 217517.110 MHz
    # CH3OCHO 48(14,35)-47(15,32)A (Methyl Formate) at 217802.057 MHz
    # CH3OH 20(1,19)-20(0,20) at 217886.39 MHz
    # CH3OCHO 45(29,16)-46(28,18)E (Methyl Formate) 217943.821 MHz
    # CH3OCHO ?? (Methyl Formate) at 218158.97 MHz
    # H2CO 3(03)-2(02)             218222.192 MHz (0.0 km/s)           yes
    # HC3N 24-23                   218324.711 MHz (-140.84 km/s)       yes
    # CH3OH 4(22)-3(12)E           218440.050 MHz (-304.79 km/s)       yes (blend with H2CO K2)
    # H2CO 3(22)-2(21)             218475.632 MHz (-348.17 km/s)       yes
    # H2CO 3(21)-2(20)             218760.071 MHz (-738.94 km/s)       yes
    # HC3N 24-23 v6=1              218854.39 MHz                       (Blend with v7=1; stronger of two vib lines?)
    # HC3N 24-23 v7=1              218860.80 MHz                       (Blend with v6=1)
    # OCS 18-17 at 218903.3555 MHz
    # HNCO 10(1,10)-9(1,9) at 218981.019 MHz
    #
    #
    # NGC4945 H2CO J=3-2 K=0,2
    #
    #my_line_list = [217.299162, 217.802057, 217.467150, 217.296605, 217.943821, 218.222192, 218.324711, 218.440050, 218.475632, 218.760071] * u.GHz
    #my_line_widths = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0] * u.km/u.s
    #my_line_names = ['CH3OH67','CH3OCHO4847','13CN21-1','13CN21-2','CH3OCHO4546','H2COJ32K0','HC3N2423','CH3OH43','H2COJ32K221','H2COJ32K210']
    # These are:
    # 13CN N=2-1 (2)               217296.605 MHz (+1732 km/s)        yes
    # CH3OH 6(15)-7(26)            217299.162 MHz (+1268.05 km/s)     yes
    # 13CN N=2-1 (1)               217467.150 MHz (+1500 km/s)        yes
    # CH3OCHO 48(14,35)-47(15,32)A 217802.057 MHz (+1045 km/s)        yes
    # CH3OCHO 45(29,16)-46(28,18)E 217943.821 MHz (+382 km/s)         yes
    # H2CO 3(03)-2(02)             218222.192 MHz (0.0 km/s)          yes
    # HC3N 24-23                   218324.711 MHz (-140.84 km/s)      yes
    # CH3OH 4(22)-3(12)E           218440.050 MHz (-304.79 km/s)      yes (blend with H2CO K2)
    # H2CO 3(22)-2(21)             218475.632 MHz (-348.17 km/s)      yes
    # H2CO 3(21)-2(20)             218760.071 MHz (-738.94 km/s)      yes
    #
    #
    # H213CO J=3-2 K=1
    #
    #my_line_list = [220.398684, 219.908486, 219.798274, 219.675114, 219.560358] * u.GHz
    #my_line_widths = [150.0, 80.0, 80.0, 80.0, 80.0] * u.km/u.s
    #my_line_names = ['13CO21','H213COJ32K1','HNCO109','HC3N2423','C18O21']
    # These are:
    # 13CO 2-1                       220398.684 MHz                     yes
    # H213CO 3(12)-2(11)             219908.486 MHz                     yes
    # HNCO 10(0,10)-9(0,9)           219798.274 MHz                     yes
    # HC3N 24-23                     219675.114 MHz                     yes
    # C180 2-1                       219560.358 MHz                     yes
    #
    #
    # H2CO J=5-4 K=1
    #
    #my_line_list = [352.199050, 351.984515, 352.603542, 351.768645, 351.6334, 351.45424, 351.236343, 351.047000, 350.905070] * u.GHz
    #my_line_widths = [40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0] * u.km/u.s
    #my_line_names = ['CH3OCHO4545','CH3OCHO3231','CH3OCHO3029','H2COJ54K1','HNCO1615','CH2NH1010','CH3OH910','Unidentified351236','CH3OH10']
    # These are:
    # CH3OCHO 45(14,32)-45(13,33)E 352199.050 MHz                     yes
    # CH3OCHO 32(1,31)-31(1,30)E   351984.515 MHz                     yes
    # CH3OCHO 30(4,27)-29(3,26)A   352603.542 MHz                     yes
    # H2CO 5(15)-4(14)             351768.645 MHz                     yes
    # HNCO 16(0,16)-15(0,15)       351633.4 MHz (+115.264 km/s)       yes (next to H2CO)              2e14
    # CH2NH 10(1,9)-10(0,10)       351454.24 MHz (+267.95 km/s)       yes (some positions)
    # CH3OH 9(5,5)-10(4,6) E       351236.343 MHz (+453.65 km/s)      yes (some positions)            1e15
    # UNIDENTIFIED                 351047.000 MHz (+615.02 km/s)      yes (strong line next to
    #                                                                      CH3OH 1-0)
    # CH3OH 1(1,1)-0(0,0) A++      350905.070 MHz (+735.98 km/s)      yes



    # Now loop over EACH line, extracting moments etc. from the appropriate region:
    # we'll also apply a transition-dependent width (my_line_widths) here because
    # these fainter lines do not have peaks as far out as the bright line.

    for line_name,line_freq,line_width in zip(my_line_names,my_line_list,my_line_widths):

        line_freq = line_freq*u.GHz
        line_width = line_width*u.km/u.s
        subcube = cube.with_spectral_unit(u.km/u.s,
                                      rest_value=line_freq,
                                      velocity_convention='optical'
                                     ).spectral_slab(peak_velocity.min()-line_width,
                                                     peak_velocity.max()+line_width)

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
                             (2*np.array(width_map)[None,:,:]**2))
        peak_sn = max_map / noisemap

        print("Peak S/N: {0}".format(np.nanmax(peak_sn)))

        # threshold at the fraction of the Gaussian corresponding to our peak s/n.
        # i.e., if the S/N=6, then the threshold will be 6-sigma
        # (this can be modified as you see fit)
        threshold = np.exp(-(peak_sn**2) / 2.)
        print("Highest Threshold: {0}".format(np.nanmax(threshold)))
        print("Lowest Threshold: {0}".format((threshold[threshold>0].min())))

        # this will compare the gaussian cube to the threshold on a (spatial)
        # pixel-by-pixel basis
        width_mask_cube = gauss_mask_cube > threshold
        print("Number of values above threshold: {0}".format(width_mask_cube.sum()))
        print("Max value in the mask cube: {0}".format(np.nanmax(gauss_mask_cube)))
        print("shapes: mask cube={0}  threshold: {1}".format(gauss_mask_cube.shape, threshold.shape))



        # this part makes a cube of velocities
        temp = subcube.spectral_axis
        velocities = np.tile(temp[:,None,None], subcube.shape[1:])

        # now we use the velocities from the brightest line to create a mask region
        # in the same velocity range but with different rest frequencies (different
        # lines)
        mask = np.abs(peak_velocity - velocities) < line_width

        # Mask on a pixel-by-pixel basis with a 3-sigma cut
        signal_mask = subcube > signal_mask_limit*noisemap

        # the mask is a cube, the spatial mask is a 2d array, but in this case
        # numpy knows how to combine them properly
        # (signal_mask is a different type, so it can't be combined with the others
        # yet - I'll add a feature request for that)
        msubcube = subcube.with_mask(mask & spatial_mask).with_mask(signal_mask).with_mask(width_mask_cube)

        # Now write output.  Note that moment0, moment1, and moment2 directories
        # must already exist...
    
        labels = {0: 'Integrated Intensity [{0}]',
                  1: '$V_{{LSR}}$ [{0}]',
                  #2: '$\sigma_v$ [{0}]',
                  2: '$FWHM$ [{0}]',
                 }

        for moment in (0,1,2):
            mom = msubcube.moment(order=moment, axis=0)
            if moment == 2:
                mom = np.multiply(2*np.sqrt(np.log(2)),np.sqrt(mom))
            hdu = mom.hdu
            hdu.header.update(cube.beam.to_header_keywords())
            hdu.header['OBJECT'] = cube.header['OBJECT']
            hdu.writeto("moment{0}/{1}_{2}_moment{0}.fits".format(moment,target,line_name), clobber=True)
            pl.figure(1).clf()
            mom.quicklook() #filename='moment{0}/{1}_{2}_moment{0}.png'.format(moment,target,line_name))
            mom.FITSFigure.colorbar.show(axis_label_text=labels[moment].format(mom.unit.to_string('latex_inline')))
            mom.FITSFigure.save(filename='moment{0}/{1}_{2}_moment{0}.png'.format(moment,target,line_name))
#
# Read parameters from dictionary
#
cubelinemoment(**params)
