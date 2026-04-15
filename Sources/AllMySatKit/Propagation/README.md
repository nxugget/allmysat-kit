# Propagation

Orbital propagation algorithms for predicting satellite positions.

## SGP4Propagator

**Complete implementation of the SGP4/SDP4 analytical orbit propagator.**

Propagates Earth-orbiting satellites from NORAD Two-Line Element (TLE) sets. Handles both near-Earth (SGP4, period < 225 min) and deep-space (SDP4, period ≥ 225 min) orbits, including geopotential resonance effects for 12-hour and synchronous orbits.

**Algorithm:**
- Secular perturbations from J₂, J₃, J₄ zonal harmonics
- Atmospheric drag via simplified general perturbation theory
- Deep-space: solar/lunar gravitational coupling (Luni-Solar perturbations)
- Resonance integration: Euler-Maclaurin numerical integrator for 12h & 24h resonant orbits
- Kepler's equation solved iteratively (Newton-Raphson, 10 iterations)
- ECI (TEME) → geodetic conversion with iterative WGS-72 oblate Earth model
- Observer look angles (azimuth, elevation, range, range-rate) via SEZ topocentric frame

**Coordinate frame:** True Equator Mean Equinox (TEME)

**References:**
- Vallado, D.A., Crawford, P., Hujsak, R., Kelso, T.S. - ["Revisiting Spacetrack Report #3"](https://celestrak.org/publications/AIAA/2006-6753/), AIAA 2006-6753
- Hoots, F.R. & Roehrich, R.L. - ["Spacetrack Report No. 3"](https://celestrak.org/NORAD/documentation/spacetrk.pdf), NORAD, 1980
- Brouwer, D. - ["Solution of the Problem of Artificial Satellite Theory Without Drag"](https://adsabs.harvard.edu/full/1959AJ.....64..378B), AJ 64, 1959

## PassPredictor

**Satellite pass prediction over a ground observer.**

Computes rise (AOS), set (LOS), and maximum elevation (TCA) times for satellite passes using a coarse-fine sweep with binary search horizon refinement.

**Algorithm:**
- Coarse scan: 30-second steps to detect horizon crossings (elevation > 0°)
- Fine scan: 1-second steps during active passes for trajectory sampling
- Horizon refinement: 15-iteration bisection on the elevation function (precision ~0.001 s)
- Pass filtering by minimum peak elevation