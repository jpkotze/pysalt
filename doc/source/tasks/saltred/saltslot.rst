.. _saltslot:

********
saltslot
********


Name
====

saltslot -- Clean raw SLOT mode data

Usage
=====

saltslot images outimages outpref gaindb xtalkfile logfile (verbose)

Parameters
==========


*images*
    String. List of images to reduce, including, if necessary, absolute or
    relative paths to the data.. Data can be provided as a comma-separated
    list, or a string with a wildcard (e.g. 'images=S20061210*.fits'), or
    a foreign file containing an ascii list of image filenames. For ascii
    list option, the filename containing the list must be provided
    preceded by a '@' character, e.g. 'images=@listoffiles.lis'. Note that
    SLOT mode fits files often contain more than one exposed frame.  In
    such cases, all frames will be reduced by default from any file
    specified in the list of images. Because the gaindb and xtalkfile
    arguments can only point to a single reference file, the images list
    must not contain data from more than one instrument. saltslot works
    specifically only on raw image data. Do not provide any files that
    have already undergone processing with SALT reduction software.

*outimages*
    String. A list of images. Data can be provided as a comma-separated
    list, or a string with a wildcard (e.g. 'outimages=rS20061210*.fits'), or
    a foreign file containing an ascii list of image filenames. For ascii
    list option, the filename containing the list must be provided
    preceded by a '@' character, e.g. 'outimages=@listoffiles.lis'. This list
    must be of the same size as the images argument list. If the
    output is intended for a different directory the absolute or relative
    path must be supplied with the file name.

*outpref*
    String. If the outpref string is non-zero in length and contains
    characters other than a blank space, it will override any value of the
    outimages argument. Output file names will use the name list provided
    in the images argument, but adding a prefix to each output file
    defined by outpref. An absolute or relative directory path can be
    included in the prefix, e.g. 'outpref=/Volumes/data/x'.

*gaindb*
    String. This is the name of an ascii table that contains the amplifier
    gains specific to a SALT instrument. The table is used to gain-correct
    amplifier raw count images. Gain is assumed to be uniform across an
    individual amplifier but assumed to vary from amplifier to amplifier.
    An example of the table format follows::

        # Database of SALTICAM CCD amplifier properties
        # 10 Aug 2006 - Telescope Data
        # READOUT GAINSTATE GAIN RDNOISE BIAS AMP
        SLOW FAINT  1.06  3.60  300 amp1
        SLOW FAINT  0.99  3.43  300 amp2
        SLOW FAINT  1.06  3.71  300 amp3
        SLOW FAINT  1.07  3.69  300 amp4
        SLOW BRIGHT 2.32  3.98  300 amp1
        SLOW BRIGHT 2.17  3.78  300 amp2
        SLOW BRIGHT 2.32  4.07  300 amp3
        SLOW BRIGHT 2.33  3.96  300 amp4
        FAST FAINT  1.55  5.22  300 amp1
        FAST FAINT  1.45  5.11  300 amp2
        FAST FAINT  1.53  5.51  300 amp3
        FAST FAINT  1.58  5.61  300 amp4
        FAST BRIGHT 4.26  6.39  300 amp1
        FAST BRIGHT 3.96  5.88  300 amp2
        FAST BRIGHT 4.21  6.35  300 amp3
        FAST BRIGHT 4.32  7.02  300 amp4

    RDNOISE refers to readout noise and BIAS refers to typical CCD bias levels.
    These data are calibrated regularly at the telescope and provided by the
    SALT project. Recent versions of the table are provided in the SALT IRAF
    distribution at salt/salticam/data/SALTICAMamps.dat and salt/pfis/data/PFISamps.dat,
    and updates will be publicized on the SALT web site at www.salt.ac.za.

*xtalkfile*
    String. This is the name of an ascii table that contains the CCD amplifier
    crosstalk coeffcients.  The table is used to subtract cross-talk contamination
    from the CCD amplifier images. Crosstalk is assumed to occur at a constant
    level across an individual amplifier, but the coefficients are assumed to
    vary across amplifer pairs. An example of the table format follows::

        # PFIS CCD amplifier crosstalk data
        # from 20041201 gain = bright distortion image, outer amps duplicated
        # Date    VCTM     2       1       4       3        6        5
        #         SRC      1       2       3       4        5        6
        2004-01-01     .001474 .001474 .001166  .001111  .001377  .001377

    The crosstalk-corrected amplifier image is given by SRC - VCTM * coeff.

*logfile*
    String. Name of an ascii file for storing log and error messages
    from the tool. The file may be new, or messages can also be appended to a
    pre-existing file.

*verbose*
    Boolean. If verbose=n, log messages will be suppressed.

Description
===========

saltslot provides basic image reduction processing for SALTICAM and
RSS SLOT mode data. SLOT mode is the continuous readout of the
detector CCDs with no deadtime overheads. It is achieved by clocking
only a small number of CCD rows towards the readout amplifiers between
each exposure.  The CCD is covered by a mask with only a thin slot cut
from the center through which the detector can expose on the sky. The
width of the slot is 144 raw pixels, which corresponds to the number
of rows clocked after each exposure. saltslot will also work on other
SALT modes and probably has value reducing Frame Transfer data.

With a maxmimum frame frequency of 10 Hz, SLOT mode is capable of
producing a large number of image frames during a track. Traditional
IRAF tools, although capable of reducing SLOT mode data, are not well
suited to large amounts of frames. IRAF is flexible and generally
robust but at the expense of economy. The IRAF approach is to perform
operations in parallel, i.e. open all files, perform one operation,
close all files, and repeat for all operations. The amount of time
performing I/O operations rapidly becomes prohibitive with large
number of images.  saltslot is a python-based script which sacrifices
sophistication and flexibility for speed. It has the ability to reduce
SLOT data faster than realtime. The algorithms used in saltslot are
generally identical to those used in other IRAF and PyRAF tools,
e.g. salt.salticam.sgain, salt.salticam.sbias, mscred.xtalkcor. Much
of the tools speed is in minimizing the amount of file I/O. Each image
is opened and closed once, and only once. All data operations are
performed in series.

The one departure from the standard IRAF algorithms is the calculation
and subtraction of the bias level. Typically, it is assumed that the
bias level varies across the detector; this assumption sacrificed for
the sake of speed. saltslot assumes that the bias level is constant
and calculates the median value in the overscan region, subtracting
this value from all image pixels. This is not usually a good
approximation, but in the case of SLOT mode, the active area of the
chip is so small that bias derivatives across the CCD rows are almost
negligible. Typically there is a one count difference between the top
row and bottom row. Investigators concerned about this difference can
remove it either by constructing and subtracting a master bias frame
or using the slower, but more robust tools salt.saltprepare,
salt.saltgain, salt.saltxtalk and salt.saltbias.

The saltslot procedure is as follows.

1. Ensure keyword and images structures are consistent with SALT
standards.

2. Gain correct images. This is simply a multiplication of all pixels by
a scalar quantity specific to amplifier, gain setting and readout
speed. The gain factor is provided by an ascii table pointed to by the
gaindb argument.

3. Subtract the bias level from images. The overscan region defined by
the BIASSEC FITS keyword is used to estimate the bias level. The median
pixel value in the overscan is adopted, and this value is subtracted
from all pixels in the image.

4. Correct for amplifier crosstalk. Each SALTICAM and RSS CCD has two
readout ammplifiers. There is crosstalk between them at the level of ~
0.1% which, provided images are not saturated or non-linear, can be
removed adequately by simple subtraction of a scaled image of one
amplifier from it's neighbour. The scaling factors are supplied as
an ascii table through the xtalkfile argument.

5. The overscan region, defined by the BIASSEC keyword, is trimmed
from each image.


Examples
========

1. To extract counts from two field stars over a series of exposures::

    --> saltslot images='S*.fits' outimages='' outpref='r'
    gaindb='/iraf/extern/salt/salticam/data/SALTICAMamps.dat'
    xtalkfile='/iraf/extern/salt/salticam/data/SALTICAMxtalk.dat'
    logfile='salt.log' verbose='yes'

Time and disk requirements
==========================

Unbinned raw slot mode images are 400KB. While saltslot processing does
result in a slightly smaller set of images after overscan stripping, it
also converts the raw images from signed 32-bit to unsigned 32-bit pixels.
This increases the size of unbinned slot mode images to 650KB, but the
advantage is a significant decrease in run time. It is recommended
to use workstations with a minimum of 512MB RAM. On a fast linux
machine with 2.8 Ghz processor and 2 Gb of RAM, one FITS file containing
two 4x4 binned slotmode exposures with a total of 1124x36 pixels can be
processed in 0.09 sec.

Bugs and limitations
====================

The current version of SALTSLOT has been tested only on SALTICAM slot
mode data. RSS SLOT mode has not yet been commissioned. Currently,
no flat fielding is performed. It is not yet clear whether
SLOT mode flat fields are useful.

Send feedback and bug reports to salthelp@saao.ac.za

See also
========

 :ref:`saltclean` :ref:`saltprepare` :ref:`saltgain` :ref:`saltbias` :ref:`saltxtalk`