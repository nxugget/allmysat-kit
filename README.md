# AllMySatKit

> **Note:** This package contains the scientific computing layer from AllMySat like orbital propagation, pass prediction, solar calculations, and coordinate transformations. For a general-purpose, maintained Swift library for SGP4/SDP4 propagation, see [SatelliteKit](https://github.com/gavineadie/SatelliteKit) by Gavin Eadie.

## What it does

- **[SGP4/SDP4 Propagation](Sources/AllMySatKit/Propagation/README.md)**: Near-earth and deep-space orbital propagation following Vallado et al. high-fidelity algorithms.
- **[TLE Parsing](Sources/AllMySatKit/TLE/README.md)**: Two-line element parsing with full validation and checksum verification.
- **[Pass Prediction](Sources/AllMySatKit/Propagation/README.md)**: Acquisition of signal (AOS) and loss of signal (LOS) with fine-tuned trajectory refinement.
- **[Sun Calculations](Sources/AllMySatKit/Astronomy/README.md)**: Sunrise, sunset, golden hour, day length, moonrise approximations using USNO methods.
- **[Visual Magnitude](Sources/AllMySatKit/Astronomy/README.md)**: Satellite brightness estimation with phase angle and Earth shadow geometry.
- **[Coordinates](Sources/AllMySatKit/Geodesy/README.md)**: DMS conversion, compass directions (8 and 16 point), grid squares per ITU-R M.431-3.
- **[Frequency Formatting](Sources/AllMySatKit/Radio/README.md)**: Radio frequency auto-scaling and band classification.

## Requirements

- Swift 5.9 or later
- iOS 14+, macOS 12+, tvOS 14+, watchOS 7+
- Foundation, CoreLocation, simd only (no external dependencies)

## Installation

Add to Package.swift as a dependency:

```swift
.package(url: "https://github.com/nxugget/allmysat-kit.git", from: "1.0.0")
```


## Quick Start

Get a satellite's position and look angles:

```swift
import AllMySatKit
import CoreLocation

let tle1 = "1 25544U 98067A   24001.50000000  .00016717  00000-0  10270-3 0  9005"
let tle2 = "2 25544  51.6416 208.5728 0006703  35.4092  75.1104 15.49560384999999"

guard let prop = SGP4Propagator(tle1: tle1, tle2: tle2),
      let eci = prop.propagate(date: Date()) else { return }

let geo = SGP4Propagator.eciToGeodetic(position: eci.position, date: Date())
print("ISS at \(geo.latitude)°, \(geo.longitude)°, \(geo.altitude) km")

let observer = CLLocation(latitude: 48.8566, longitude: 2.3522)
let look = SGP4Propagator.lookAngles(
    position: eci.position, velocity: eci.velocity,
    observer: observer, date: Date()
)
print("Azimuth \(look.azimuth)°, elevation \(look.elevation)°")
```

Predict passes over the next 7 days:

```swift
let observer = CLLocation(latitude: 48.8566, longitude: 2.3522)
let passes = PassPredictor.predictPasses(
    propagator: prop,
    observer: observer,
    from: Date(),
    to: Date(timeIntervalSinceNow: 7 * 86400),
    minElevation: 10.0
)

for pass in passes {
    print("TCA at \(pass.tca.time), max elevation \(pass.maxElevation)°")
}
```

Get sunrise and golden hour times:

```swift
let sun = SunCalculator.compute(
    latitude: 48.8566,
    longitude: 2.3522,
    date: Date()
)

print("Sunrise: \(sun.sunrise)")
print("Golden hour: \(sun.goldenHourStart)")
print("Day length: \(sun.dayLengthHours) hours")
```

## Testing

Run all 130 tests across 10 suites:

```bash
swift test
```

Covers SGP4 propagation, TLE parsing, solar calculations, coordinate systems, grid squares, frequency formatting, and performance stress tests.

## License

Apache License 2.0
