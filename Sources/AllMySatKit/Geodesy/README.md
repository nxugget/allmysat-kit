# Geodesy

Coordinate systems and geodesic calculations for Earth.

## CoordinateConverter

**Coordinate transformation utilities for satellite tracking.**

Converts between decimal degrees, degrees-minutes-seconds (DMS), and 3D direction vectors. Provides compass direction quantization and geodesic distance conversions.

**Functions:**
- Decimal ↔ DMS conversion
- Azimuth/elevation → unit direction vector (SIMD3) for 3D scene rendering
- 16-point and 8-point compass rose quantization
- Distance unit conversion (km ↔ mi ↔ nm)

**References:**
- DMS standard: [ISO 6709](https://en.wikipedia.org/wiki/ISO_6709) - "Standard representation of geographic point location by coordinates"
- Nautical mile: 1 nm = 1852 m (International Hydrographic Organization definition)
- Compass rose: 32-point system, here using 16 and 8-point subsets

## GridSquareCalculator

**Maidenhead Grid Square Locator System (QTH Locator).**

Implements the grid square system used worldwide by amateur radio operators and satellite ground stations. Converts between geographic coordinates and Maidenhead locator strings at four precision levels.

**Encoding scheme (ITU-R M.431-3):**
- Field (2 chars, A–R): 20° longitude × 10° latitude
- Square (2 digits, 0–9): 2° × 1°
- Subsquare (2 chars, a–x): 5′ × 2.5′
- Extended square (2 digits, 0–9): 30″ × 15″

**Also provides:**
- Great-circle distance (via CoreLocation's Vincenty implementation)
- Initial bearing between locators (forward azimuth on the WGS-84 ellipsoid)
- Grid square bounds for map overlays
- Locator validation (regex against ITU specification)

**References:**
- [ISO 6709](https://en.wikipedia.org/wiki/ISO_6709) - Geographic coordinate representation
- Maidenhead grid system used in amateur radio and satellite operations
