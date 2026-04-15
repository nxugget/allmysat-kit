# Astronomy

Solar position and satellite visibility calculations.

## SunCalculator

**Solar position and sunrise/sunset computation.**

Calculates sunrise, sunset, and day length for any location and date using the Earth's axial tilt and the equation of time. The lunar phase is approximated using the synodic month period.

**Algorithm:**
- Solar declination: δ = −23.44° × cos(360°/365 × (d + 10))
- Hour angle from zenith angle: cos(ω) = (sin(−0.833°) − sin(φ)sin(δ)) / (cos(φ)cos(δ))
- Solar noon: 12:00 − longitude/15°
- Lunar cycle: 29.53-day synodic month approximation

**References:**
- NOAA Solar Calculator - [gml.noaa.gov/grad/solcalc](https://gml.noaa.gov/grad/solcalc/)
- Meeus, J. - ["Astronomical Algorithms"](https://dn710207.ca.archive.org/0/items/astronomicalalgorithmsjeanmeeus1991/Astronomical%20Algorithms-%20Jean%20Meeus%20%281991%29.pdf), 2nd ed., Willmann-Bell, 1998, Chapters 12–15

## VisualMagnitudeCalculator

**Apparent visual magnitude and naked-eye visibility for satellites.**

Determines whether a satellite pass is visible to the naked eye by computing apparent brightness at sampled points along the trajectory, accounting for observer darkness, satellite illumination, range, and phase angle.

**Algorithm:**
- Sun elevation: NOAA-based J2000.0 model - mean anomaly, equation of center, ecliptic longitude, right ascension, hour angle via GMST, then elevation from the fundamental astronomical triangle
- Earth shadow: cylindrical shadow model using Earth angular radius as seen from satellite altitude, with altitude correction
- Apparent magnitude: m = m_std + 5 log₁₀(r/r_std) + Δm_phase
  - Range correction: inverse square law (5 mag per decade of distance)
  - Phase correction: Lambertian diffuse sphere model - Δm = −2.5 log₁₀[(1+cos θ)/(1+cos θ_std)]
- Visibility criterion: apparent magnitude ≤ 6.0 (naked-eye limit in dark skies)

**References:**
- McCants, M. - ["Visual Satellite Observer's FAQ: Predicting Satellite Brightness"](https://www.heavens-above.com/faq.aspx)
- Astronomical magnitude system: Pogson, N.R. - ["Magnitudes of Thirty-six of the Minor Planets"](https://academic.oup.com/mnras/article/17/1/12/956950), MNRAS 17, 1857
- Standard magnitude convention: 1000 km range, 50° phase angle
