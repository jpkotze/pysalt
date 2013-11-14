#! /usr/bin/python
# -*- coding: utf-8 -*-
#
# Author                   Version      Date
# -----------------------------------------------
# Keith Smith (Nottingham)   1.0        19 July 2011

# Module for finding SALT guide stars
# - full changelog is in /doc/changelog.txt
# - manual is (will be?) in /doc/manual.txt

global __version__
__version__ = "1.0"
__date__ = "19 July 2011"
__author__ = "Keith Smith"
__copyright__ = "Copyright "+ __author__ + " " + __date__

__doc__="\nSALT guide star finder, version "+__version__ +"""

Finds guide stars suitable for use with the Southern African Large Telescope
Uses the HST Guide Star Catalogue 2.3.2 and the VizieR catalogue service

Usage: python guide_stars.py [OPTIONS] [TARGET]

TARGET should be the coordinates of the target or a SIMBAD-resolvable name
Acceptable formats are:
	colon-seperated sexagesimal eg. 15:43:17.2 -18:56:12.9
	space-seperated sexagesimal eg. '15 43 17.2' '-18 56 12.9'
	decimal hours and degrees eg. 15.7837 -18.9543
	target name eg. m31 or 'RV Cen'
All coordinates should be J2000 right ascension and declination

OPTIONS are as follows, arguments are compulsory for both long and short forms:
		--help              Prints this help
	-v	--verbose           Prints useful information during execution
	-d	--debug             Prints debugging information, implies -v
	-c	--current           Uses the current telescope pointing
	-o	--ocs				Use the current OCS targets, writes the result to OCS
								Don't try to combine this with -c, -t or command line input
	-t	--target-id=ID      Queries the database for the specified target ID,
								writes the result to the database
	-f	--filter=FILTER     Uses the FILTER filter,
		                        defaults to R
	-i	--instrument=INS    Uses the INS instrument,
		                        defaults to RSS
	-r	--radius=RADIUS     Assumes target size is RADIUS arcseconds in radius,
		                        defaults to 2 arcseconds"""


# import required modules

import urllib2              # reading URLs
import StringUtil           # converting decimal/sexagesimal
from xml.dom import minidom # XML parsing
import csv                  # CSV parsing
import StringIO
from numpy import *         # arrays
import sys
import getopt               # command line switches
import MySQLdb              # SQL queries

# define exceptions

class GuideStarError(Exception):
    pass    # misc errors in this module
class VizError(Exception):
    pass # errors thrown by VizieR
class BadInput(Exception):
    StringUtil
    pass # incorrect input

# define global variables, later read from the command line

verbose = False
debug = False

# define helper functions

def usage():
    print __doc__
    raise SystemExit(2)		# 2 is the UNIX code for bad command line input apparently

def isDecimal(string):
    "returns true if the input string converts to a valid decimal"
    try:
        float(string)
        return True
    except ValueError:
        return False

def vprint(string):
    "prints the string if verbose mode is on"
    global verbose
    if verbose == True:
        print string

def dprint(string):
    "prints the string if debug mode is on"
    global debug
    if debug == True:
        print string

def checkInput(targetRA, targetDec, filter, instrument, targetRadius):
    "validates input to the main function"

    name_mode = False		# unless we find otherwise below
    # removing this line breaks the unit test, due to direct call to this function

    # parse RA and dec into sexagesimal or target name for passing to VizieR

    print targetRA, targetDec
    if targetDec == '':
        # in target name mode
        name_mode = True
        
    elif isDecimal(targetDec):
        # dec is decimal
        if isDecimal(targetRA)==False:  # check RA in same format
            raise BadInput, 'RA and dec appear to be in different formats'  
        targetDec = StringUtil.dmsStrFromDeg(float(targetDec))
        targetRA = StringUtil.dmsStrFromDeg(float(targetRA))
        print targetRA, targetDec
    elif isDecimal(targetRA):
        # RA is in decimal, but dec wasn't
        # pretty sure there are no valid target names which are decimals...
        # need to check this as decimals pass as valid sexa!
            raise BadInput, 'RA and dec appear to be in different formats'
 
    elif StringUtil.checkDMSStr(targetDec) == True:
        # target already colon-seperated sexa
        # couldn't check this first as this function returns True
        # for decimals (for some reason), but False if space seperated sexa
        if StringUtil.checkDMSStr(targetRA) == False:   # check RA in same format
            raise BadInput, 'RA and dec appear to be in different formats'
            pass    # both already in correct format
    elif StringUtil.checkDMSStr(targetDec.replace(' ',':')) == True:
        # is a valid space seperated sexa, convert to colon seperated
        if (StringUtil.checkDMSStr(targetRA.replace(' ',':')) == False or
            (targetRA.replace(' ',':')==targetRA)):
        # check RA works colon seperated and wasn't to start with
            raise BadInput, 'RA and dec appear to be in different formats'

        targetDec=targetDec.replace(' ',':')
        targetRA=targetRA.replace(' ',':')
    else:
        # dec isn't decimal, or space/colon seperated sexa, or blank
        raise BadInput, 'Format of declination not recognised, was given: '+str(targetDec)

    # having now got our ra and dec into a common format, check they are sensible numbers

    if name_mode==False:

    # Splits the RA into a list with several elements - works for both sexa and dec
        checkRA = StringUtil.splitDMSStr(targetRA)
        print checkRA
        if checkRA[0]=='-':	# RA is -ve
            raise BadInput, 'RA is negative'
        elif len(checkRA)==3:	# decimal
            if float(targetRA)>360 or float(targetRA)<0: # RA decimal isn't the range 0-360 degrees
                raise BadInput, 'RA appears to be in decimal, but isn\'t in the range 0 to 360 degrees'
        elif len(checkRA)==5: #sexa
            
            if (int(checkRA[1])>23 or int(checkRA[1])<0 # RA hours aren't in the range 0-23
            or int(checkRA[2])>59 or int(checkRA[2])<0 # RA minutes arent in the range 0-59
            or int(checkRA[3])>59 or int(checkRA[3])<0): # RA seconds arent in the range 0-59
                raise BadInput, 'RA appears to be in sexagesimal, but isn\'t in the range 0:0:0 to 23:59:59'
        else: # not sexa or decimal, but survived above
            raise GuideStarError, 'RA passed valid formatting checks, but not in sexa or decimal. This shouldn\'t happen.'

        # Same checks on Dec
        checkDec = StringUtil.splitDMSStr(targetDec)
        if len(checkDec)==3:	# decimal
            if abs(float(targetDec))>90: # Dec decimal isn't the range -90-90 degrees
                raise BadInput, 'Dec appears to be in decimal, but isn\'t in the range -90 to +90 degrees'
        elif len(checkDec)==5:	#sexa
            if (abs(int(checkDec[1])>90) # Dec degrees aren't in the range -90-90
            or int(checkDec[2])>59 or int(checkDec[2])<0 # Dec arcminutes arent in the range 0-59
            or int(checkDec[3])>59 or int(checkDec[3])<0): # Dec arcseconds arent in the range 0-59
                raise BadInput, 'RA appears to be in sexagesimal, but isn\'t in the range 0:0:0 to 23:59:59'
        else: # not sexa or decimal, but survived above
            raise GuideStarError, 'RA passed valid formatting checks, but not in sexa or decimal. This shouldn\'t happen.'

        if name_mode==True:
            vprint('Using target name ' + targetRA)
        else:
            vprint('Target coordinates are ' + targetRA + ' ' + targetDec)

    # parse declination to deal with idiosyncrasies of VizieR and urllib2

    if name_mode==False:
        if "+" in targetDec:    # source is northern
            targetDec = "%2b" + targetDec
        elif "-" in targetDec:  # source is southern
            targetDec = "%20" + targetDec
        else:                   # no sign, assume +ve and add it
            targetDec = "%2b+" + targetDec

    else:
        targetRA = targetRA.replace(' ','%20')
        # spaces break the URL for some reason

    # convert filters to those in the GSC

    if filter == "":
        vprint('No filter specified, defaulting to photographic R_F')
        filter = 'F'
    elif filter in ('j', 'F', 'V', 'N'):
        vprint("Selected photographic " + filter + " filter")
        pass    # these are the native filters in the GSC
    elif filter in ('U', "u'", 'u'):
        vprint('Specified '+ filter + ' filter')
        filter = 'U'	# there IS a U field in GSC, but hardly any entries
    elif filter in ('B', "g'", 'b', 'v'):
        vprint('Specified ' + filter + ' filter, closest GSC band is photographic Bj')
        filter = 'j'   # there is also a Bmag field in GSC, but less entries
    elif filter in 'y':
        vprint('Specified ' + filter + ' filter, closest GSC band is photographic V')
        filter = 'V'
    elif filter in ('R', "r'"):
        vprint('Specified ' + filter + ' filter, closest GSC band is photographic F')
        filter = 'F'
    elif filter in ('I', "i'", "z'"):
        vprint('Specified ' + filter + ' filter, closest GSC band is photographic N')
        filter = 'N'
    else:
        raise BadInput, "Filter '%s' not recognised" % str(filter)

    # check instrument input

    instrument = instrument.lower()
    if instrument == "":
        vprint('No instrument specified, defaulting to RSS')
        instrument = 'rss'
    elif instrument in ['rss','pfis']:
        vprint('Selected RSS')
        instrument = 'rss'
    elif instrument in ["salticam","scam"]:
        vprint('Selected SALTICAM; SALTICAM values not yet validated')
        instrument = 'scam'
    elif instrument == "hrs":
        vprint('Selected HRS; HRS values not yet validated')
        instrument = 'hrs'
    else:
        raise BadInput, 'Instrument "' + str(instrument) + '" not recognised'

    # check radius

    if targetRadius == '':
        vprint('No target radius specified, defaulting to 2 arcsec')
        targetRadius = 2.
    elif 0. < targetRadius <= 240.:
        vprint('Selected target radius of ' + str(targetRadius) + ' arcsec')
    elif targetRadius > 240.:
        raise BadInput, 'Target radius '+str(targetRadius)+' arcsec is larger than the science FoV'
    else:
        raise BadInput, 'Target radius of ' + str(targetRadius) + 'arcsec is invalid'

    return (targetRA, targetDec, filter, instrument, targetRadius)

def queryUrl(url):
    "accesses the input url and returns the response"

    request = urllib2.Request(url)
    opener = urllib2.build_opener()
    global __version__	# this doesn't seem to work, not sure why. Not vital.
    request.add_header('User-Agent', 'SALT guide star finder/' +  __version__  + ' www.salt.ac.za')
    dprint('User-Agent: SALT guide star finder/' + __version__ + ' www.salt.ac.za')
    try:
        data = opener.open(request).read()
    except:
        raise GuideStarError, "Could not connect to VizieR. Please check your internet connection"
    return data

def queryTelescopePointing():
    "Queries the current telescope pointing, returns current RA and dec"

    tcs_url = "http://icd.salt/xml/salt-tcs-icd.xml"

    vprint("Querying current telescope pointing")
    dprint("Querying URL " + tcs_url)

    tcs_response = queryUrl(tcs_url)
    tcs_xml = minidom.parseString(tcs_response)

    # The block below uses some hard coded childNode indices. This is because each element has a \n as its first element. These will break if the XML gets reformatted...

	# Look in the at the <Cluster> elements
    for cluster in tcs_xml.getElementsByTagName("Cluster"):
        # Look for the current pointing cluster
        if cluster.childNodes[1].firstChild.data == "tcs pointing actual info":
            # Inside those are <DBL> elements
            for dbl in cluster.getElementsByTagName("DBL"):
                if dbl.childNodes[1].firstChild.data == "RA":
                    ra = dbl.childNodes[3].firstChild.data # note in radians
                if dbl.childNodes[1].firstChild.data == "Dec":
                    dec = dbl.childNodes[3].firstChild.data # note in radians	

	# TCS returns RA and Dec in radians. Convert to decimal degrees
    ra = 180.0 * float(ra) / pi
    dec = 180.0 * float(dec) / pi

    # Make sure RA is in range 0-360
    if ra < 0.:
        ra=ra+360.

	# Convert RA from decimal degrees to decimal hours
	ra = 24.*ra/360.

	vprint("Current pointing is RA = %.3f, Dec = %.3f" % (ra, dec))
	# No need to force these into a nice format, as the checkInput function will do that for us later

	return ra, dec


def constructVizUrl(targetRA, targetDec, min_r, max_r, filter, min_mag, max_mag):
    "constructs the url for the VizieR guide star query"

    url = "http://vizier.u-strasbg.fr/cgi-bin/asu-xml?" # base URL for VizieR XML queries
    url += "-source=I/305/out"                  # use HST GSC 2.3.2
    url += "&-c=" + targetRA
# old dec line
#    url += "%20" + targetDec + "&-c.eq=J2000"   # target section
    url += targetDec + "&-c.eq=J2000"   # target section
    url += "&-c.rm=%s,%s" % (min_r, max_r)      # anulus range (arcmin)
    url += "&-out=%smag&%smag=%s..%s" % (filter, filter, min_mag, max_mag)
    url += "&-out.max=1000"                        # max entries
    url += "&-out=_r,_RA*-c.eq,_DE*-c.eq"       # calculate 2000 RA, dec and distance
    url += ",Class&Class=0"                     # only objects flagged as stars
    url += "&-oc.form=sexa"                     # output coords in sexagesimal
    url += "&-sort=-_r"                         # sort by decreasing r
    url += "&-mime=CSV"                         # data in CSV format

    return url



def parseVizResponse(response):
    "parses the response from VizieR into into a list of data for each star"

# extract CSV table from XML

    if "****" in response:
        # Between 2007 and 2011 Vizier changed how it reported errors.
        # This is reinstated kludgey code to find them, probably a better way to do this
        errordata=""
        for line in response:
            if "****" in response:
                errordata += line	# actually this appends the WHOLE response...
        raise VizError, "Vizier generated an error. The error returned was:\n" + errordata


    xml = minidom.parseString(response)

    table = xml.getElementsByTagName('CSV')
    if len(table)==0:
        info = xml.getElementsByTagName('INFO')
        for element in info:
            id = element.getAttribute('ID')

            # Note this is for the way Vizier reported errors in 2007. By 2011, this had changed (see above and changelog)
            if id=='Errors':
                #found an error
                raise VizError, 'Vizier generated an error. The error was:' + element.firstChild.data
        # no CSV table was found ie. no data
        # needs some error checking - might be a bad url or target
        return [],0,[]
    colsep = table[0].getAttribute('colsep')        # CSV columns seperator
    headlines = table[0].getAttribute('headlines')  # CSV header lines
    table = table[0].firstChild.data                # the CSV table itself

    # convert to from unicode
    colsep = colsep.encode('ASCII')
    headlines = int(headlines)


# extract data from CSV table

    csvreader=csv.reader(StringIO.StringIO(table), delimiter=colsep)


# this next bit is very kludgy, weird csvreader object and array manipulations
    for row in csvreader:   # can't index csvreader!
        if len(row)>0:  # blank rows in the reader!
            if 'csvtable' in locals():      # checks to see if csvtable exists
            # python has no exist() function!
                csvtable=vstack((csvtable,array([row])))
            # fragile, dimension mismatch possible
            else:
                csvtable=array([row])
            # can't append to blank arrays!

    n_stars=len(csvtable)-headlines

# split off the Class column, as an added bonus won't break if there isn't one
    index = (csvtable[0,:]!='Class')  # boolean index of first row != Class
    csvtable = csvtable[:,index]        # just retain columns in the index

    vprint("Retrieved " + str(n_stars )+ " records from VizieR")

    return csvtable, n_stars, headlines



def sortResults(table, n_stars, headlines, min_r, max_r, min_mag, max_mag):
    "Sorts the results to select the best guide stars"

    vprint('Sorting...')

    # split the header off the table
    header=table[:headlines]
    data=table[headlines:]

    top=header[0,:]
    dprint('Columns are: ' + str(top))

    r_index = top=='_r'
    if not True in r_index:
        raise GuideStarError, ('Could not find a radius column in search results\n' +
		  'Returned columns are: ' + str(top))
    dprint('Identified radius column ' + str(r_index))
    r=data[:,r_index]

    mag_index = (top=='Vmag') + (top=='Nmag') + (top=='Bmag') +\
			    (top=='Fmag') + (top=='Umag') + (top=='jmag')
    # hard coded brute force
    if not True in mag_index:
        raise GuideStarError, ('Could not find a magnitude column in search results\n' +
		  'Returned columns are: ' + str(top))
    dprint('Identified magnitude column ' + str(mag_index))
    mag=data[:,mag_index]

    # convert to floats
    r=r.astype(float64)
    mag=mag.astype(float64)

    r_pref=select([r<min_r, (min_r<=r) & (r<=max_r), r>max_r],
			      [min_r-r, 0., r-max_r])

    mag_pref=select([mag<min_mag, (min_mag<=mag) & (mag<=max_mag), mag>max_r],
			      [min_mag-mag, 0., mag-max_mag])

    pref = 2*r_pref + mag_pref	# weighted such that 1 mag = 0.5 arcmin

    sort_index = lexsort(keys=(mag.T,pref.T),axis=1)	# decide the order from the pref list
	# in order, so LOWEST first
    # note the index is crazy - lists the indexes in order, not the order of indices!
    sorted = data[sort_index]
    pref_sort = pref[sort_index]

    n_stars = min(n_stars,6)	# return at most 6 stars

    # there must be an easier way to print these strings!
    dprint('Preferences calculated, best: '+str(pref_sort[0,0].tolist()[0])+
	       ', worst: ' +str(pref_sort[0,-1].tolist()[0])+
	       ', cutoff: '+str(pref_sort[0,n_stars-1].tolist()[0])+
	       ', returning '+str(n_stars)+' stars')

    # replace some of the column headings
    top[top.tolist().index('_r')]='Offset'
    top[top.tolist().index('_RAJ2000')]='RA (J2000)'
    top[top.tolist().index('_DEJ2000')]='Dec (J2000)'

    header[0,:]=top

    # now passes the header separately
    ## put the header back on
    #data = vstack((header,sorted[:,:n_stars][0]))
    data = sorted[:,:n_stars][0]

    return data, header, n_stars


# main function

def findGuideStars(targetRA, targetDec="",
					filter='', instrument='', targetRadius=2.):

# doc string

    '''Finds SALT-suitable guide stars based upon the HST GSC 2.3.2.
    Pass target RA and Dec (J2000), and optionally
    filter, instrument, and target radius in arcsec.
    All input should be as strings.
    Keith Smith, Nov 2007'''

    status = 'Preparing...'
    vprint(status)

    (targetRA, targetDec, filter, instrument, targetRadius)= \
	    checkInput(targetRA, targetDec, filter, instrument, targetRadius)

# set up instrument-specific settings
# this should go into the check input function, but too much passing
# for the moment - will change later
# alternatively use classes?
    instrument = instrument.lower()
    if instrument == 'rss':
        pref_max_r = 4.9     # prefered outer edge of anulus in arcmin
        pref_min_r = 4.     # prefered inner edge of anulus in arcmin
        abs_max_r = 5.      # absolute limit on outer edge
        abs_min_r = 1.      # absolute limit on inner edge
        pref_max_mag = 10.  # prefered magnitude limits
        pref_min_mag = 9.
        abs_max_mag = 19.5   # absolute magnitude limits
        abs_min_mag = 7.
    elif instrument == 'scam':
        # salticam guider has not yet been installed - use values from rss
        # if you fill this in, please change warning message in input check above
        pref_max_r = 4.9     # prefered outer edge of anulus in arcmin
        pref_min_r = 4.     # prefered inner edge of anulus in arcmin
        abs_max_r = 5.      # absolute limit on outer edge
        abs_min_r = 1.      # absolute limit on inner edge
        pref_max_mag = 10.  # prefered magnitude limits
        pref_min_mag = 9.
        abs_max_mag = 19.5   # absolute magnitude limits
        abs_min_mag = 7.
    elif instrument == "hrs":
        # hrs has not yet been built - use values from rss
        # if you fill this in, please change warning message in input check above
        pref_max_r = 4.9     # prefered outer edge of anulus in arcmin
        pref_min_r = 4.     # prefered inner edge of anulus in arcmin
        abs_max_r = 5.      # absolute limit on outer edge
        abs_min_r = 1.      # absolute limit on inner edge
        pref_max_mag = 10.  # prefered magnitude limits
        pref_min_mag = 9.
        abs_max_mag = 19.5   # absolute magnitude limits
        abs_min_mag = 7.
    else:
        raise GuideStarError, "Instrument settings not found, but passed input checking. This shouldn't happen."


# modify radii if the source is very large

    targetRadius = (targetRadius / 60.) # convert from arcsec to arcmin
    if targetRadius > pref_max_r:
        if targetRadius > abs_max_r:
            return "Target is larger than the guide star field of view"
        else:
            pref_max_r = abs_max_r
    if targetRadius > abs_min_r:
        abs_min_r = targetRadius    # + a bit to avoid vignetting?
        if abs_min_r > pref_min_r:
            pref_min_r = abs_min_r


# set up query
    url = constructVizUrl(targetRA, targetDec, abs_min_r, abs_max_r, filter, abs_min_mag, abs_max_mag)

    vprint("Querying...")
    dprint("Using URL:\n" + url)

    # retrieve data
    response = queryUrl(url)
    data, n_stars, header = parseVizResponse(response)

    if n_stars>0:
        data, header, n_stars = sortResults(data, n_stars, header,
							pref_min_r, pref_max_r, pref_min_mag, pref_max_mag)


    if n_stars>1:
        vprint('Selected best ' + str(n_stars) + ' guide stars')
        status = 'Found ' + str(n_stars) + ' stars'
    elif n_stars==1:
        vprint('Found one guide star')
        status = 'Found 1 star'
    elif n_stars==0:
        status = 'Failed, no suitable stars found'


    return n_stars, data, header

# end main function


def commandLine(argv):
    # executes if module is run from the command line

    dprint("Reading command line options")

    # read command line options
    try:
        opts,args = getopt.getopt(sys.argv[1:],"vdct:f:i:r:o",
                ["verbose","debug","current", "target-id=","filter=","instrument=","radius=","ocs","help"])
    except getopt.GetoptError, inst:
        print inst
        print 'Use --help to get a list of options'
        sys.exit(2)

    ra, dec, filter, ins, radius, target_id = "","","","","",""
    use_current_pointing = False
    use_ocs = False
    global verbose
    global debug

    # parse them to the relevant variables
    for opt, arg in opts:
        if opt in ('--help'):
            usage()
        elif opt in ('-v','--verbose'):
            verbose=True
        elif opt in ('-d','--debug'):
            verbose=True	# implied
            debug=True
        elif opt in ('-f','--filter'):
            filter = arg
        elif opt in ('-i','--instrument'):
            ins = arg
        elif opt in ('-r','--radius'):
            radius = float(arg)
        elif opt in ('-t','--target-id'):
            target_id = arg
        elif opt in ('-c','--current'):
            use_current_pointing = True
        elif opt in ('-o','--ocs'):
            use_ocs = True
        else:
            print 'Unknown option: ' + opt
            usage()

    for argument in args:
        if ra=="":
            ra = argument
        elif dec=="":
            dec = argument
        else:
            #too many arguments
            raise BadInput, 'Too many arguments, takes one or two input arguments'

    if use_current_pointing == True:
        ra,dec = queryTelescopePointing()

    if use_ocs == True:
        vprint("Reading from current OCS table")

        # Nasty MySQL querying
        ocs_db=MySQLdb.connect(host="sdbdev", user="guide_stars", passwd="guide_stars4ocs", db="sdb_v5a")
        cursor=ocs_db.cursor()
        cursor.execute("select OcsTargetInformation_Id, RaRad, DecRad, GuideStarRa, GuideStarDec from OcsTargetInformation where GuideStarRa=0 && GuideStarDec=0")
        ocs_table=cursor.fetchall()
        cursor.close()

        dprint("Response from SQL query was:")
        dprint(ocs_table)

        if len(ocs_table)==0:
            vprint("Nothing to change")

        for line in ocs_table:
            (ocs_id,ra_rad,dec_rad,gs_ra_rad,gs_dec_rad)=line
            # Convert from radians to DMS
            ra=StringUtil.dmsStrFromDeg((ra_rad*180./pi)/24.) # extra /24 is to give RA in hours
            dec=StringUtil.dmsStrFromDeg(dec_rad*180./pi)
            try:
                n_stars,data,header = findGuideStars(ra,dec,filter,ins,radius)

            except:
                vprint("Error in target id %f, trying next target" % ocs_id)
                continue	# skip this target

            dprint(data)

            # Select the 'best' guide star. Kludgy
            gs_ra_hms=data[0,2]
            gs_dec_hms=data[0,3]
            vprint("For OCS_Id %i, guide star RA %s Dec %s" % (ocs_id, gs_ra_hms, gs_dec_hms))

            # Convert to radians
            gs_ra_deg=StringUtil.degFromDMSStr(gs_ra_hms.replace(' ',':'))*24. # *24 because was in hms, not dms
            gs_dec_deg=StringUtil.degFromDMSStr(gs_dec_hms.replace(' ',':'))

            gs_ra_rad=gs_ra_deg*pi/180.
            gs_dec_rad=gs_dec_deg*pi/180.

            # Write those results to the OCS database using mySQL
            cursor=ocs_db.cursor()	# remembers ocs_db from earlier, so no need to keep setting
            SQL_string="update OcsTargetInformation set GuideStarRa=%s, GuideStarDec=%s where OcsTargetInformation_Id = %s" % (gs_ra_rad,gs_dec_rad,ocs_id)
            dprint("SQL string is %s" % SQL_string)
            cursor.execute(SQL_string)
            ocs_db.commit()	# have to do this to actually run the query!
            cursor.close()

        return "","",""
        #ocs_array = array(ocs_table)
        #
        #ocs_ids=ocs_array[:,0]
        #ra_rads=ocs_array[:,1]
        #dec_rads=ocs_array[:,2]
        ## Guide star RA and Dec are empty, or query would not have returned

    if target_id != "":
        raise BadInput, "Not implemented yet"
        ra,dec = queryDatabase(target_id)

    if ra=="":	# no target was specified
        raise BadInput, 'No target specified'

    n_stars,data,header = findGuideStars(ra,dec,filter,ins,radius)

    dprint(n_stars)

    # Fall back to other filters if no stars found
    # Order is based on GSC coverage in the Southern hemisphere
    # Down to 19.5 mag, coverage order is (R, I, B, V, U), see 2008AJ....136..735L
    while n_stars==0:
	# U is the last filter, so if it failed this all hope is lost (was pretty bad if we even tried U)
        if filter == 'U':
            vprint('Nothing else to fall back to. Try modifying your query')
            status = 'Failed, no suitable stars found'
            break

    # Cycling through the filters, remember closest bands are I->N, B->j
        else:
            vprint("Failed, found nothing in " + filter + ", trying other filters")
            if filter in ('F',""):	# blank defaults to R
                filter = 'N'
            elif filter == 'N':
                filter = 'j'
            elif filter == 'j':
                filter = 'V'
            elif filter == 'V':
                filter = 'U'
            vprint("Falling back to " + filter + " filter...")
            n_stars, data, header = findGuideStars(ra, dec, filter, ins, radius)



    return n_stars,data,header

# module input

if __name__ == "__main__":

    # pass stuff to the command line parser
    status,data,header = commandLine(sys.argv)
    dprint(status)
    vprint(header)
    print data
    sys.exit(0)
