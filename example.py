from pysolar.solar import *
import datetime
from datetime import datetime as dt
from datetime import timedelta
from pysolar.util2 import get_sunrise_sunset_transit

########################## Sunrise Sunset ##################################################

SUN_ZENITH_SUNRISE_SUNSET = 0 # Sunrise sunset
SUN_ZENITH_CIVIL_TWILIGHT = 1 # Civil twilights
SUN_ZENITH_NAUTICAL_TWILIGHT = 2 # Nautical twilights
SUN_ZENITH_ASTRONOMICAL_TWILIGHT = 3 # Astronomical twilights

#ADD LONGITUDE-LATITUDE
obsLon = -17.778  #lon is negative to the west.
obsLat = 28.606

date = dt.now()
tz = pytz.timezone('UTC')
date = date.replace(tzinfo=tz)

#SUNRISE SUNSET

sr, ss, tr = get_sunrise_sunset_transit(obsLat, obsLon, date, SUN_ZENITH_SUNRISE_SUNSET)
print("ASTRONOMICAL TWILIGHT")
print('sunrise: ', sr.strftime("%Y-%m-%d %H:%M:%S"))
print('sunset:', ss.strftime("%Y-%m-%d %H:%M:%S"))
print('transit:', tr.strftime("%Y-%m-%d %H:%M:%S"))


#ASTRONOMICAL TWILIGHT
sr, ss, tr = get_sunrise_sunset_transit(obsLat, obsLon, date, SUN_ZENITH_ASTRONOMICAL_TWILIGHT)

print("ASTRONOMICAL TWILIGHT")
print('sunrise: ', sr.strftime("%Y-%m-%d %H:%M:%S"))
print('sunset:', ss.strftime("%Y-%m-%d %H:%M:%S"))
#print('transit:', tr.strftime("%Y-%m-%d %H:%M:%S"))