![license](https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square)

# SunMoonCalculator
A very simple Sun/Moon calculator without using JPARSEC library.

This is a Python port of the excellent [ephemerides (in Java)](http://conga.oan.es/~alonso/doku.php?id=blog:sun_moon_position) by Tomás Alonso Albi.

For specified observer date/time and location (Geodetic: longitude/latitude), following information is calculated:
- **Sun**
  - Rise, set, noon (transit) time.
  - Azimuth, elevation and transit elevation angles.
  - Distance.
- **Moon**
  - Rise, set, transit time.
  - Azimuth, elevation and transit elevation angles.
  - Distance.
  - Age.
  - Moon illumination (percentage).
  - Moon phase.
  
