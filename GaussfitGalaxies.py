"""
A tool for fitting Gaussians to fit 2D Gaussians to sources specified
in a catalog/region file.
"""
import regions
from astropy import units as u
#import paths
from gaussfit_catalog import gaussfit_catalog, gaussfit_image
from astropy.table import Table,Column

def tryint(x):
    try:
        return int(x)
    except:
        return x

def data_to_table(fit_data):
    names = fit_data.keys()
    numnames = [ii for ii,nm in enumerate(names)]
    stripnames = [nm for nm in names]
    stripnames = [fullname for nnm,fullname in sorted(zip(numnames,stripnames))]
    names = [fullname for nnm,fullname in sorted(zip(numnames,names))]
    namecol = Column(name='Name', data=stripnames)
    colnames = ['amplitude', 'center_x', 'center_y', 'fwhm_major', 'fwhm_minor', 'pa',
                'deconv_fwhm_major', 'deconv_fwhm_minor', 'deconv_pa',
                'chi2', 'chi2/n', 'e_amplitude', 'e_center_x', 'e_center_y',
                'e_fwhm_major', 'e_fwhm_minor', 'e_pa', 'ampguess', 'peak',
                'success',]
    columns = [Column(name=k, data=[fit_data[entry][k].value
                                    if hasattr(fit_data[entry][k],'value')
                                    else fit_data[entry][k]
                                    for entry in names],
                      unit=(fit_data[names[0]][k].unit
                            if hasattr(fit_data[names[0]][k], 'unit')
                            else None))
               for k in colnames]

    tbl = Table([namecol]+columns)

    tbl['ampguess'].description = 'Amplitude guess (peak intensity - median background'
    tbl['peak'].description = 'Peak intensity'

    return tbl



if __name__ == "__main__":

    for regfn, contfn, name in (
        ('Meier2015NGC253Positions-Updated.reg', 'NGC253-H2COJ32K02_H2COJ32K0_moment0_widthscale1.0_sncut2.0_widthcutscale1.0.fits', 'NGC253-H2COJ32K02-H2COJ32K0GaussFit'),
        #('Meier2015NGC253Positions-Updated.reg', 'NGC253-H2COJ32K02_H2COJ32K210_moment0_widthscale1.0_sncut2.0_widthcutscale1.0.fits', 'NGC253-H2COJ32K02-H2COJ32K210GaussFit'),
        #('Meier2015NGC253Positions-Updated.reg', 'NGC253-H2COJ54K1_H2COJ54K1_moment0_widthscale1.0_sncut2.0_widthcutscale1.0.fits', 'NGC253-H2COJ54K1-H2COJ54K1GaussFit'),
        #('Meier2015NGC253Positions-Updated.reg', 'NGC253-H2COJ54K23_H2COJ54K321_moment0_widthscale1.0_sncut2.0_widthcutscale1.0.fits', 'NGC253-H2COJ54K23-H2COJ54K321GaussFit'),
    ):

        regs = regions.read_ds9(regfn)

#        contfn = 'gaussfit/'+contfn

        fit_data = gaussfit_catalog(contfn, regs, radius=1.0*u.arcsec,
                                    prefix=name+"_",
                                    max_radius_in_beams=5,
                                    max_offset_in_beams=2,
                                    #raise_for_failure=True,
                                    savepath='gaussfit')

        tbl = data_to_table(fit_data)

        tbl.rename_column("chi2/n", "chi2_n")
        tbl.write("gaussian_fit_table_{0}.ipac".format(name),
                  format='ascii.ipac',
                  overwrite=True)
