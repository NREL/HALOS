
"""
HALOS Solar Positioning Module 

Solar Position - Azimuth and Zenith 

"""

import csv
import datetime

import numpy

import pysolar.solar 

"""
The script reads the weather file and returns Solar Azimuth angle, Solar Zenith Angle for a complete year.
The script uses pysolar -- https://pysolar.readthedocs.io/en/latest/

The results obtained can be compared with the results available at 
https://www.esrl.noaa.gov/gmd/grad/solcalc/azel.html
"""

def reading_weather_file(weather_file):
        weather_data = {"dni": numpy.zeros(8760, dtype=float), 
                        "year": numpy.zeros(8760, dtype=int), 
                        "month": numpy.zeros(8760, dtype=int), 
                        "day": numpy.zeros(8760, dtype=int), 
                        "hour": numpy.zeros(8760, dtype=int)}
      
        
        reader = csv.reader(open(weather_file, 'r'))
        # First line is column headers for categorical data

        header_keys = next(reader)
        # Second line is data for first line of headers; record lat and lon
        fline = next(reader)
        d = {}
        for i in range(len(fline)):
            d[header_keys[i]] = fline[i]
            
        weather_data["lat"] = float(d["Latitude"])
        weather_data["lon"] = float(d["Longitude"])
        weather_data["time_zone"] = int(d["Time Zone"])
        
        # Third line is column headers for time-series data, just track dni
        keys = next(reader)
        year_idx = keys.index("Year")
        month_idx = keys.index("Month")
        day_idx = keys.index("Day")
        hour_idx = keys.index("Hour")
        if d["Source"] == "IWEC":
            dni_idx = keys.index("Beam")
        else:
            dni_idx = keys.index("DNI")
        #Record time-series data
        t = 0
        for line in reader:
            weather_data["dni"][t] = float(line[dni_idx])
            weather_data["year"][t] = int(line[year_idx])
            weather_data["month"][t] = int(line[month_idx])
            weather_data["day"][t] = int(line[day_idx])
            weather_data["hour"][t] = int(line[hour_idx])
            
            t += 1
            if t == 8760:
                break
        return weather_data

def get_single_solar_position(weather_data,idx):
    """
    Get Solar Position for single period

    Parameters
    ----------
    weather_data : Dataframe
    idx : hour id 

    Returns
    -------
    zenith : Zenith Angle of Sun
    azimuth : Azimuth Angle of Sun

    """
    lon = - weather_data["lon"]
    lat = weather_data["lat"]
    year = weather_data["year"][idx]
    month = weather_data["month"][idx]
    day = weather_data["day"][idx]
    hour = weather_data["hour"][idx]
    time_zone = weather_data["time_zone"]
    time_of_the_year = datetime.datetime(year,month,day,hour,tzinfo=datetime.timezone.utc) + datetime.timedelta(hours=time_zone)
    zenith = float(90) - pysolar.solar.get_altitude(lat, lon, time_of_the_year)     #Zenith is the difference between 90 deg and elevation
    azimuth = pysolar.solar.get_azimuth(lat, lon, time_of_the_year)
    return zenith, azimuth

def solar_position(weather_data):
    """
    Get solar position for complete year. 

    Parameters
    ----------
    weather_data : Dataframe

    Returns
    -------
    Zenith : Zenith of every hour
    Azimuth : Azimuth of every hour
    """
    Zenith = []
    Azimuth = []
    lon = - weather_data["lon"]
    lat = weather_data["lat"]
    time_zone = weather_data["time_zone"]
    for idx in range(len(weather_data["year"])):
        year = weather_data["year"][idx]
        month = weather_data["month"][idx]
        day = weather_data["day"][idx]
        hour = weather_data["hour"][idx]
        time_of_the_year = datetime.datetime(year,month,day,hour,tzinfo=datetime.timezone.utc) + datetime.timedelta(hours=time_zone)
        zenith = float(90) - pysolar.solar.get_altitude(lat, lon, time_of_the_year)     #Zenith is the difference between 90 deg and elevation
        Zenith.append(zenith)
        azimuth = pysolar.solar.get_azimuth(lat, lon, time_of_the_year)
        Azimuth.append(azimuth)
        
    return (Zenith, Azimuth)

if __name__ == "__main__":
    import flux_model
    weather_data = flux_model.ReadWeatherFile("./../weather_files/USA NV Tonopah Airport (TMY3).csv")
    (Zenith, Azimuth) = solar_position(weather_data)
    print(Azimuth)
    print("Zenith is")
    print(Zenith)
    