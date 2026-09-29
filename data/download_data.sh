#! /bin/sh

#
# Download the relevant data files on
#  - tilt angle,
#  - modulation potential from NM data,
#  - smoothed sunspot number, and
#  - heliospheric magnetic field amplitude and solar wind speed (27-day averages) from OMNIWeb.
#
# Existing data files will be backed up before overwriting.
#

set -e

# end date for OMNIWeb query
OMNIWEB_ENDDATE="20260916"

tilt_angle()
{
    echo '====== Tilt angle ======'
    if [ -f angle.dat ]
    then
	cp angle.dat angle.bak
    fi

    wget -N http://wso.stanford.edu/Tilts.html
    cat > angle.dat <<EOF
# WSO Computed "Tilt Angle" of the Heliospheric Current Sheet
# http://wso.stanford.edu/Tilts.html
#
# Carr Rot    Start Date	R_av	R_n	R_s	L_av	L_n	L_s
EOF
    # extract only the data lines
    awk '/^CR/' Tilts.html >> angle.dat
}

neutron_monitor()
{
    echo '====== NM data ======'
    if [ -f nm.dat ]
    then
	cp nm.dat nm.bak
    fi

    wget -N https://cosmicrays.oulu.fi/phi/Phi_mon.txt
    cat > nm.dat <<EOF
# Modulation parameter (in MV) reconstructed from ground-based CR data
# https://cosmicrays.oulu.fi/phi/Phi_mon.txt
#
# Year Jan   Feb   Mar   Apr   May   Jun   Jul   Aug   Sep   Oct   Nov   Dec  Annual
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
EOF
    # extract only the data lines, remove a trailing whitespace
    awk '/^[0-9][0-9 ]+/' Phi_mon.txt | sed 's/ \r$//;' >> nm.dat
}

sunspot_number()
{
    echo '====== SSN ======'
    if [ -f ssn.dat ]
    then
	cp ssn.dat ssn.bak
    fi

    wget -N https://www.sidc.be/SILSO/DATA/SN_ms_tot_V2.0.txt
    cat > ssn.dat <<EOF
# 13-month smoothed monthly total sunspot number
# https://www.sidc.be/SILSO/DATA/SN_ms_tot_V2.0.txt
# Column 1-2: Gregorian calendar date (year, month)
# Column 3: Date in fraction of year for the middle of the corresponding month
# Column 4: Monthly smoothed total sunspot number.
# Column 5: Monthly mean standard deviation of the input sunspot numbers.
# Column 6: Number of observations used to compute the corresponding monthly mean total sunspot number.
# Column 7: Definitive/provisional marker. A blank indicates that the value is definitive. A '*' symbol indicates that the monthly value is still provisional and is subject to a possible revision (Usually the last 3 to 6 months)
#
EOF
    cat SN_ms_tot_V2.0.txt >> ssn.dat
}

omniweb()
{
    echo '====== OMNIweb data ======'
    if [ -f omniweb.dat ]
    then
	cp omniweb.dat omniweb.bak
    fi

    wget --post-data  "activity=retrieve&res=27day&spacecraft=omni2_27day&start_date=19640101&end_date=${OMNIWEB_ENDDATE}&vars=8&vars=24&scale=Linear&ymin=&ymax=&view=0&charsize=&xstyle=0&ystyle=0&symbol=0&symsize=&linestyle=solid&table=0&imagex=640&imagey=480&color=&back=" https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi -O omniweb_download.txt
    # extract only the data lines
    awk '/^[0-9][0-9\. ]+/' omniweb_download.txt > omniweb.dat
}

tilt_angle
neutron_monitor
sunspot_number
omniweb
