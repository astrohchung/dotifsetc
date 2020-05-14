# DOTIFS Exposure Time Calculator (S/N Calculator)

## NAME: 
    dotifsetc()
  
## PURPOSE:
  DOTIFS exposure time calculator (ETC, or S/N calculator) is developed to 
  provide an expected signal to noise ratio information in various observing
  conditions and targets to community users who are observing or planing
  to do observation with Devasthal Optical Telescope Integral Field
  Spectrograph (DOTIFS).

## REQUIRED PACKAGES:
  packages with lower than below versions would also work.
  
    PYTHON: 3.6.5 
    NUMPY: 1.14.3
    SCIPY: 1.1.0
    MATPLOTLIB: 2.2.2
    DOTIFSETC_UTIL: included in the distribution.

## INSTALLATION:
  type and execute the following commands on the git-installed Linux 
  terminal
  
    git clone git@github.com:astrohchung/dotifsetc.git
    
  if you would like to check whether the program is successfully installed,
  please run the simple example code as below.
  
    cd dotifsetc
    python3 dotifsetc_example_simple.py
  
  above code will produce two result files: 'dotifs_etc.ps' and '3600.ps'. 
  the results should be identical to 'dotifs_etc_ex.ps' and '3600_ex.ps' in
  this distribution.

## CALLING SEQUENCE:
    from dotifsetc import dotifsetc 
    result=dotifsetc(exptime=3600, band='r', magnitude=20., skymagnitude=22,
                     oname='dotifs_snr.ps', source='obj_sc', z=0., stype=0, 
                     wstep=(3700./3000), pixel=None, ltrghost=False, 
                     soc=False, skysub=True, wavearr=None, inputflux=None, 
                     inputwave=None, show=False, save=False, 
                     plotrange=[3700, 7400], run=True)
    
    print(result.wave)
    print(result.sourceflux)
    print(result.snr)
    print(result.signal)
    result.show=True
    result.plot()
    
  dotifsetc() without any parameters will also work and give a result with
  predefined observing condition.
    
## INPUT PARAMETERS:
    EXPTIME: exposure time in seconds (default: 900 seconds)
    BAND: photometric band which will be used to calculate source flux. 
        Only SDSS ugriz bands are supported. (default: 'r')
    MAGNITUDE: source magnitude at the defined photometric band. The template 
        flux will be scaled based on this input magnitude. (AB magnitude)
        (default: 17 mag)
    SKYMAGNITUDE: sky magnitude at the defined photometric band. The template 
        flux will be scaled based on this input magnitude. (AB magnitude)
        (default: 17 mag)
    ONAME: name of output file (default: 'dotifs_snr.ps')
    SOURCE: target name. ETC read SED of selected target from target
        templates or generate SED based on target option. source parameter is
        comprised of two strings - source type and value. available source 
        type and value is listed as below. (default: 'obj_sc')
        source types(source values):
            obj(s0, sa, sb, sc, bulge, elliptical, starb1, starb2, starb3,
                starb4, starb5, starb6): model spectrum of various targets.
                data is obtained from kc96.
                (http://www.stsci.edu/hst/observatory/crds/cdbs_kc96.html)
            arcflat(Xenon): spectrum of Xenon arc lamp. data from Newport
                catalog.
            wavecal(Kr, HgNe, KrHgNe): spectrum of wavelength calibration
                sources.
            sky: spectrum of sky at various geometry between Sun, Moon, Earth 
                and target. source value of this type of object is provided 
                separately using stype parameter.
            blackbody(temparature in Kelvin): blackbody spectrum which follows 
                Planck distribution.
            const: flux with constant magnitude.
            Example: 
                sa galaxy: dotifsetc(source='obj_sa')
                elliptical galaxy: dotifsetc(source='obj_elliptical')
                arcflat: dotifsetc(source='arcflat_Xenon')
                sky: dotifsetc(source='sky', stype=0)
                     dotifsetc(source='sky', stype='p0_p135_n90_p45')
                blackbody: dotifsetc(source='blackbody_5500')
                const: dotifsetc(source='const')
    Z: redshift of the source. target SED will be redshifted according to this
        input parameter. It works only when source type is obj.
        (default: 0.000)
    STYPE: sky type selecting parameter. 
        format: mssep_mtsep_malt_talt
            the format indicates spectrum of sky at various geometry
            between Sun, Moon, Earth and target. all in degrees. data is
            obtained from ESO sky calculator.
            (https://www.eso.org/observing/etc/skycalc/) Predefined options 
            are listed in the description of STYPE parameter. 
            mssep(moon-sun separation), mtsep(moon-target separation),
            malt(moon altitude), talt(target altitude)
            there are more options on model sky spectrum, and users may
            use one with their preferred option by putting model sky file
            (in photon count) at sky_spectrum directory and modify 
            sky_templates.fmt file to let etc to understand the file.
        user can choose model either by option index or model name. 
        currently, only four sky template spectrum models are supported as a 
        below list. (default: 0)
        model options:
            p0_p135_n90_p45: newmoon, target altitude = 45 degrees.
            p90_p90_p45_p45: halfmoon. 
                moon and target altitude = 45 degrees.
                separation between moon and target = 90 degrees.
            p180_p90_p45_p45: fullmooon. 
                moon and target altitude = 45 degrees.
                separation between moon and target = 90 degrees.
    AIRMASS: airmass of the target. This airmass value is used to calculate
        the sky transmission using the sky extinction data (default: 1)
    WSTEP: wavelength step size of the output spectra in Angstrom.
        (default: 1.233 angstrom)
    PIXEL: wavelength step size in pixel unit. if this pixel parameter is 
        defined, then wstep parameter will be ignored.
        (default: None)
    LTRGHOST: set this keyword to on/off littrow ghost on top of the source 
        spectrum. (default: False)
    SOC: set this keyword to on/off second order contamion on top of the 
        source spectrum. (default: False)
    SKYSUB: set this keyword to show skysubtracted result. (default: False)
        when this keyword is True, output S/N, signal, and noise count is
        sky-subtracted result. user can obtain identical result by manually
        calculate the result from non-skysubtracted observation result and
        sky observation result. 
        Example:
            wosky=dotifsetc(source='obj_sc', magnitude=20, skymagnitude=22, 
                            skysub=True) #sky subtracted
            wsky=dotifsetc(source='obj_sc', magnitude=20, skymagnitude=22,
                           skysub=False) #sky non-subtracted result
            sky=dotifsetc(source='sky', magnitude=22, skysub=False)
            wosky.signal=wsky.signal-sky.signal
            wosky.noise=(wsky.noise**2+sky.noise**2)**0.5
    WAVEARR: user can use their own wavelength grid by providing wavelength
        vecton in numpy array format. WSTEP and PIXEL parameters will not be
        used when this parameter is provided. (default: None)
    INPUTFLUX: user can use user defined source flux as a source of ETC 
        in unit of erg/s/cm^2/Ang. inputwave should be provided as well.
        (default: None)
    INPUTWAVE: wavelength of inputflux in angstrom. (default: None)
    SKYTRANS: set this keyword to the bandname of the sky transmission file
         as written in the response_curves.fmt file, if a user wants to use
         sky transmission file instead of sky extinction file. 
         For example, if set this keyword as 'default', then the sky 
         transmission is calculated from the extinction data 
         ('sky_extinction.dat'). The data format of the extinction data 
         should be 'wavelength' at column 1 and 'mag/airmass' at column 2.
         If set this keyword as 'skytrans', then the 
         'sky_transmission_eso_sky_calc.dat' will be used. AIRMASS keyword
         will be ignored in this case.
         (default: False)
    RMEDSN: set this keyword to return an median SN value among band. 
        Range of band is determined by BAND parameter. Two bands are 
        supported: 'g' and 'r' (default: False)
    TSCALE: user can scale the entire throughput of the optics by changing
        this parameter (default: 1)
    SHOW: set this keyword to view the result with matplotlib window.
        (default: False)
    SAVE: set this keyword to save the plot in the output file. 
        (default: False)
    PLOTRANGE: set the wavelength range of the result plot in two elements
        list format. (wavelength in angstrom) (default: [3700,7400])
    RUN: set this keyword to calculate the result when dotifsetc class is 
        generated. (default: True)

## ADDITIONAL INPUT PARAMETERS:
    user can modify below parameters as attributes of the dotifsetc class.
    user should execute dotifsetc.run() function to get the new result.
    .ITPKIND: choose interpolation method to interpolate the parameter vector
        on the output wavelength grid. available options are listed as below.
        (default: 'cubic') 
        options: 'linear', 'nearest', 'zero','slinear', 'quadratic', 'cubic', 
        'previous', 'next'
        detail of each method is described in the scipy documentation. 
        (scipy.interpolate.interp1d)
    .PRI: telescope primary mirror diameter in meter (default: 3.6)
    .SEC: telescope secondary mirror diameter in meter (default: 0.915)
    .SOURCESAMPLINGSIZE: size of one spatial element in square arcsecond. 
        Useful to set this as 1 if one wants to use this exposure time
        calculator for single object. For exmaple, the calculator can be
        used to estimate the S/N of MOS target, by setting this keyword as
        1, MAGNITUDE as a fiber magnitude and SKYSAMPLINGSIZE as a field of
        view of fiber.
        (default: 0.831384) (size of hexagon with 0.4 arcsecond side)
    .SKYSAMPLINGSIZE: set this keyword to the size of spatial element. 
        if None, the value will be the same as SOURCESAMPLINGSIZE.
        (default: None) 
    .DISPERSION: dispersion of the spectrograph in Angstrom per micron.
        (default: 0.082222) 
    .PIXELSIZE: size of CCD pixel in micron. (default: 15)
    .NPIX_SPA: size of PSF on CCD along spatial direction in pixel unit. 
        this number is used to calculate readout noise count at each 
        wavelength bin. (default: 5)
    .RN: readout noise count in ADU (default: 2)
    .DARK: dark current in unit of ADU per an hour per pixel. (default: 0)
    
## OUTPUT PARAMETERS:
    user can read output parameters as attributes of the dotifsetc class.
    .SOURCEFLUX: SED of input template in a unit of erg/cm2/sec/angstrom.
    .WAVE: wavelength of SED
    .SNR: signal to noise ratio
    .SIGNAL: electron signal count of the given observation.
    .NOISE: electron noise count of the given observation. this is
        quadrature summation of noises from source, sky, readout, and
        dark current. (when skysub=True, sky noise is added one more 
        time.)

## EXAMPLES:
    from dotifsetc import dotifsetc
    to view the result of below examples, user should save them in a
    variable as result=dotifsetc(), or set show=True, or set save=True 
    - Test run:
        dotifsetc(show=True, save=True)
        matplotlib window will be poped up, and the result will be saved
        in 'dotifs_snr.ps' file.
    - r band surface brightness=17, Sc type galaxy at redshift z=0.04, 
    10 minutes exposure with output name of 'Sc_rmag17_z0.04.ps':
        dotifsetc(exptime=600, band='r', magnitude=17, source='obj_sc', 
                  z=0.04, oname='Sc_rmag17_z0.04.ps', save=True)
    - Sa type galaxy, apply second order contamination, and Littrow 
    ghost. change wavelength bin size as 2.5 pixels:
        dotifsetc(source='obj_sa', soc=True, ltrghost=True, pixel=2.5)
    - Wavelength Calibration source spectrum. (Krypton lamp):
        dotifsetc(source='wavecal_Kr', exptime=5)
    - Arc lamp source spectrum. (Xenon lamp):
        dotifsetc(source='arcflat_Xenon', exptime=5)
    - sa type galaxy, plotrange from 5000-6800 angstrom:
        dotifsetc(source='obj_sa', plotrange=[5000,6800])
    - sa type galaxy, g band surface birightness=18, with fullmoon sky
        and skymagnitude at g band=20
        dotifsetc(source='obj_sa', band='g', magnitude=18, stype=2, 
                  skymagnitude=20)
    - sa type galaxy, read output data.
        result=dotifsetc(source='obj_sc')
        snr=result.snr
        signal=result.signal
        noise=result.noise
        wavelength=result.wave
    - run ETC after modify some optional parameters. (primary diameter
        of the telescope as 8 meters, change interpolation scheme as 
        linear)
        result=dotifsetc(run=False)
        result.pri=8
        result.itpkind='linear'
        result.save=True
        result.run()

## NOTE:
    the signal count does not account for contribution from noise. 
    (noise is not added)
    
## MODIFICATION HISTORY:
    v1.0.0: Haeun Chung, 2015, Jun. 28, IUCAA, First version
    v1.1.0: Haeun Chung, 2016, Jun. 9, IUCAA, Modifed for internal 
            distribution
    v1.2.0: Haeun Chung, 2017, Oct. 24, SNU, Add calibration sources. 
            Add Arc lamp source. Fix usd_asahi=1. Add Asachi filter at
            AOI=10 deg.
    v1.3.0: Haeun Chung, 2018, Nov. 12, SNU, Add Littrow ghost. Change 
            gen_cal scheme to add non-zero flux to every pixel
    v2.0.0: Haeun Chung, 2018. Aug. 8, Translated from IDL to Python 
            and tested against original version.
    v2.1.0: Haeun Chung, 2020. May 14, Steward Observatory (Day 57 of 
            WFH due to COVID-19)
            - add airmass parameter. Now sky transmission can be 
            calculated from sky extinction data (mag/airmass). Direct 
            use of sky transmission data is also available. 
            - add rmedsn parameter
            - add tscale parameter
            - add sourcesampling and skysampling parameter
            - change interpolation method of sky template file. Now 
            log10(skyflam) value is interpolated to avoid
            negative output from cubic interpolation.