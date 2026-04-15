# TLE

NORAD Two-Line Element set parsing and orbital parameter extraction.

## TLEParser

**Parser for the NORAD/USSPACECOM Two-Line Element set format.**

Extracts epoch, classical orbital elements, and derived Keplerian parameters from standard TLE strings. The TLE format encodes satellite state in a fixed-width 69-character-per-line text format defined by NORAD.

**Parsed fields (line 2):**
- Inclination (columns 9–16)
- Right Ascension of Ascending Node (columns 18–25)
- Eccentricity (columns 27–33, implicit leading "0.")
- Argument of Perigee (columns 35–42)
- Mean Anomaly (columns 44–51)
- Mean Motion (columns 53–63, revolutions/day)

**Epoch parsing (line 1, columns 19–32):**
- Format: YYDDD.DDDDDDDD
- Year disambiguation: 00–56 → 2000s, 57–99 → 1900s

**Derived parameters (Kepler's Third Law):**
- Semi-major axis: a = (μ/n²)^(1/3), where n is mean motion in rad/s
- Orbital period: T = 1440/n_rev minutes
- Apogee altitude: a(1+e) − R_earth
- Perigee altitude: a(1−e) − R_earth

**References:**
- Kelso, T.S. - ["CelesTrak: NORAD Two-Line Element Set Format"](https://celestrak.org/columns/v04n03/)
- Hoots, F.R. & Roehrich, R.L. - ["Spacetrack Report No. 3"](https://celestrak.org/NORAD/documentation/spacetrk.pdf)
