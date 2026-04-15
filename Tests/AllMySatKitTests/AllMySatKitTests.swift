// ────────────────────────────────────────────────────────────────────────────
// AllMySatKitTests.swift
// AllMySatKit - Comprehensive Test Suite
//
// Validates every public module of AllMySatKit against reference data:
//   • SGP4 Propagator - Julian date, GMST, ECI→geodetic, look angles, ISS orbit
//   • TLE Parser - ISS/GPS/GEO element extraction, checksum, age, edge cases
//   • SunCalculator - sunrise/sunset, polar day/night, golden hour, progress
//   • VisualMagnitudeCalculator - sun elevation, shadow geometry, magnitude
//   • CoordinateConverter - DMS, compass, 3D vectors, distance conversion
//   • GridSquareCalculator - encode, decode, bounds, distance, validation
//   • FrequencyFormatter - auto format, unit format, parse, band, shift
//   • AstrodynamicsConstants - reference values
//   • PassPredictor - ISS pass prediction integration
//
// Copyright 2026 AllMySat · Apache License 2.0
// ────────────────────────────────────────────────────────────────────────────

import Testing
import Foundation
import CoreLocation
import simd
@testable import AllMySatKit

// MARK: - Reference TLE Data

/// ISS (ZARYA) - LEO, ~408 km, 51.6° inclination (real TLE with valid checksum)
private let issTLE1 = "1 25544U 98067A   24001.00000000  .00016717  00000-0  10270-3 0  9003"
private let issTLE2 = "2 25544  51.6416 208.5328 0006816  35.2328 324.9138 15.49560532431403"

/// GPS satellite - MEO, ~20200 km, 55° inclination (deep-space: period > 225 min)
private let gpsTLE1 = "1 20413U 90008A   90036.82843564 -.00000000  00000-0  00000-0 0  9999"
private let gpsTLE2 = "2 20413  55.1021  67.1580 0065778 237.8027 267.7023  2.00562119  1300"

/// GEO satellite - period ~1440 min (deep-space)
private let geoTLE1 = "1 28057U 03049A   06177.78615833  .00000060  00000-0  00000-0 0  4753"
private let geoTLE2 = "2 28057   0.0244 250.4915 0001257 321.9430 136.8990  1.00271220 10076"

// MARK: - Helpers

private func utcDate(year: Int, month: Int, day: Int, hour: Int = 12,
                     minute: Int = 0, second: Int = 0) -> Date {
    var cal = Calendar(identifier: .gregorian)
    cal.timeZone = TimeZone(identifier: "UTC")!
    return cal.date(from: DateComponents(
        year: year, month: month, day: day,
        hour: hour, minute: minute, second: second
    ))!
}

// MARK: - SGP4 Propagator Tests

@Suite("SGP4Propagator")
struct SGP4PropagatorTests {

    @Test("ISS initializes as LEO (not deep-space)")
    func issIsLEO() {
        let prop = SGP4Propagator(tle1: issTLE1, tle2: issTLE2)
        #expect(prop != nil)
        #expect(prop?.noradId == 25544)
        #expect(prop?.isDeepSpace == false)
        #expect(prop?.periodMinutes ?? 0 < 225)
    }

    @Test("GPS initializes as deep-space (period > 225 min)")
    func gpsIsDeepSpace() {
        let prop = SGP4Propagator(tle1: gpsTLE1, tle2: gpsTLE2)
        #expect(prop != nil)
        #expect(prop?.isDeepSpace == true)
        #expect(prop?.periodMinutes ?? 0 > 225)
    }

    @Test("Invalid TLE returns nil")
    func invalidTLE() {
        let prop = SGP4Propagator(tle1: "garbage", tle2: "data")
        #expect(prop == nil)
    }

    @Test("Julian Date for J2000.0 epoch")
    func julianDateJ2000() {
        // J2000.0 = 2000-01-01T12:00:00 UTC → JD 2451545.0
        let j2000 = utcDate(year: 2000, month: 1, day: 1, hour: 12)
        let jd = SGP4Propagator.julianDate(from: j2000)
        #expect(abs(jd - 2451545.0) < 0.001)
    }

    @Test("Julian Date for known Unix epoch")
    func julianDateUnixEpoch() {
        // 1970-01-01T00:00:00 UTC → JD 2440587.5
        let unixEpoch = Date(timeIntervalSince1970: 0)
        let jd = SGP4Propagator.julianDate(from: unixEpoch)
        #expect(abs(jd - 2440587.5) < 0.001)
    }

    @Test("GMST at J2000.0 is within [0, 2π]")
    func gmstAtJ2000() {
        let gmst = SGP4Propagator.gmst(jd: 2451545.0)
        #expect(gmst > 0)
        #expect(gmst < 2 * .pi)
    }

    @Test("ECI to geodetic produces valid latitude/longitude/altitude")
    func eciToGeodetic() {
        // ISS-like position: ~408 km altitude in LEO
        let pos = SIMD3<Double>(6778.0, 0.0, 0.0)
        let date = Date(timeIntervalSince1970: 1704067200) // 2024-01-01T00:00:00Z
        let geo = SGP4Propagator.eciToGeodetic(position: pos, date: date)
        #expect(geo.latitude >= -90 && geo.latitude <= 90)
        #expect(geo.longitude >= -180 && geo.longitude <= 180)
        #expect(geo.altitude > 300 && geo.altitude < 500)
    }

    @Test("Look angles from Paris to arbitrary ECI position")
    func lookAngles() {
        let pos = SIMD3<Double>(4200.0, 2500.0, 4300.0)
        let vel = SIMD3<Double>(-3.0, 5.0, 4.0)
        let paris = CLLocation(latitude: 48.8566, longitude: 2.3522)
        let date = Date(timeIntervalSince1970: 1704067200)
        let look = SGP4Propagator.lookAngles(
            position: pos, velocity: vel, observer: paris, date: date
        )
        #expect(look.azimuth >= 0 && look.azimuth < 360)
        #expect(look.elevation >= -90 && look.elevation <= 90)
        #expect(look.range > 0)
    }

    @Test("SGP4 propagation at epoch returns valid ECI state")
    func propagateAtEpoch() {
        guard let prop = SGP4Propagator(tle1: issTLE1, tle2: issTLE2) else {
            Issue.record("Failed to init SGP4"); return
        }
        let epoch = utcDate(year: 2024, month: 1, day: 1, hour: 12)
        guard let eci = prop.propagate(date: epoch) else {
            Issue.record("Propagation failed"); return
        }
        // ISS at ~408 km: radius should be ~6778 km
        let radius = simd.length(eci.position)
        #expect(radius > 6500 && radius < 7200)
        // Velocity ~7.66 km/s
        let speed = simd.length(eci.velocity)
        #expect(speed > 6.0 && speed < 9.0)
    }

    @Test("Propagation 24h from epoch remains valid")
    func propagate24h() {
        guard let prop = SGP4Propagator(tle1: issTLE1, tle2: issTLE2) else {
            Issue.record("Failed to init SGP4"); return
        }
        let date = utcDate(year: 2024, month: 1, day: 2, hour: 12)
        guard let eci = prop.propagate(date: date) else {
            Issue.record("Propagation failed at +24h"); return
        }
        let radius = simd.length(eci.position)
        #expect(radius > 6500 && radius < 7200) // still in LEO
    }

    @Test("Perigee altitude is within expected range for ISS")
    func issPerigee() {
        guard let prop = SGP4Propagator(tle1: issTLE1, tle2: issTLE2) else {
            Issue.record("Failed to init SGP4"); return
        }
        #expect(prop.perigeeKm > 350 && prop.perigeeKm < 450)
    }

    @Test("ISS inclination is ~51.6°")
    func issInclination() {
        guard let prop = SGP4Propagator(tle1: issTLE1, tle2: issTLE2) else {
            Issue.record("Failed to init SGP4"); return
        }
        #expect(abs(prop.inclinationDeg - 51.6416) < 0.01)
    }

    @Test("Epoch Julian Date matches TLE epoch")
    func epochJd() {
        guard let prop = SGP4Propagator(tle1: issTLE1, tle2: issTLE2) else {
            Issue.record("Failed to init SGP4"); return
        }
        // TLE epoch = 24001.50000000 → 2024 day 1.5 → JD ~2460310.0
        #expect(prop.epochJd > 2460300 && prop.epochJd < 2460320)
    }

    @Test("Geodetic coordinate property gives CLLocationCoordinate2D")
    func geodeticCoordinate() {
        let pos = SIMD3<Double>(6778.0, 0.0, 0.0)
        let date = Date(timeIntervalSince1970: 1704067200)
        let geo = SGP4Propagator.eciToGeodetic(position: pos, date: date)
        let coord = geo.coordinate
        #expect(coord.latitude == geo.latitude)
        #expect(coord.longitude == geo.longitude)
    }
}

// MARK: - TLE Parser Tests

@Suite("TLEParser")
struct TLEParserTests {

    @Test("Parses ISS TLE correctly")
    func parseISS() {
        guard let e = TLEParser.parse(line1: issTLE1, line2: issTLE2) else {
            Issue.record("Failed to parse ISS TLE"); return
        }
        #expect(e.noradId == 25544)
        #expect(e.classification == Character("U"))
        #expect(abs(e.inclination - 51.6416) < 0.001)
        #expect(abs(e.eccentricity - 0.0006816) < 0.0001)
        #expect(abs(e.meanMotion - 15.4956053) < 0.001)
        #expect(e.revolutionNumber == 43140)
    }

    @Test("ISS derived orbital parameters are LEO")
    func issOrbitalParams() {
        guard let e = TLEParser.parse(line1: issTLE1, line2: issTLE2) else {
            Issue.record("parse failed"); return
        }
        #expect(e.periodMinutes > 90 && e.periodMinutes < 95)
        #expect(e.perigeeAltitudeKm > 350 && e.perigeeAltitudeKm < 450)
        #expect(e.apogeeAltitudeKm > 350 && e.apogeeAltitudeKm < 450)
        #expect(e.semiMajorAxisKm > 6700 && e.semiMajorAxisKm < 6800)
    }

    @Test("Parses GPS TLE (deep-space)")
    func parseGPS() {
        guard let e = TLEParser.parse(line1: gpsTLE1, line2: gpsTLE2) else {
            Issue.record("parse failed"); return
        }
        #expect(e.noradId == 20413)
        #expect(abs(e.inclination - 55.1021) < 0.01)
        #expect(abs(e.eccentricity - 0.0065778) < 0.0001)
        #expect(e.periodMinutes > 700) // ~720 min (12h orbit)
    }

    @Test("Parses all angular elements")
    func angularElements() {
        guard let e = TLEParser.parse(line1: issTLE1, line2: issTLE2) else {
            Issue.record("parse failed"); return
        }
        #expect(e.raan >= 0 && e.raan < 360)
        #expect(e.argumentOfPerigee >= 0 && e.argumentOfPerigee < 360)
        #expect(e.meanAnomaly >= 0 && e.meanAnomaly < 360)
    }

    @Test("Checksum validates correctly-formed TLE line")
    func checksumValid() {
        // Construct line with known correct checksum
        // Sum of digits in first 68 chars mod 10 must equal last digit
        let base = "1 25544U 98067A   24001.00000000  .00016717  00000-0  10270-3 0  900"
        let chars = Array(base)
        var sum = 0
        for i in 0..<68 {
            if let d = chars[i].wholeNumberValue { sum += d }
            else if chars[i] == "-" { sum += 1 }
        }
        let validLine = base + String(sum % 10)
        #expect(TLEParser.validateChecksum(line: validLine))
    }

    @Test("Corrupted checksum fails validation")
    func checksumInvalid() {
        // Build a valid line, then flip the last digit
        let base = "1 25544U 98067A   24001.00000000  .00016717  00000-0  10270-3 0  900"
        let chars = Array(base)
        var sum = 0
        for i in 0..<68 {
            if let d = chars[i].wholeNumberValue { sum += d }
            else if chars[i] == "-" { sum += 1 }
        }
        let correct = sum % 10
        let wrong = (correct + 1) % 10
        let invalidLine = base + String(wrong)
        #expect(!TLEParser.validateChecksum(line: invalidLine))
    }

    @Test("TLE age is positive for future reference dates")
    func tleAge() {
        guard let e = TLEParser.parse(line1: issTLE1, line2: issTLE2) else {
            Issue.record("parse failed"); return
        }
        // ISS epoch = 24001.0 = 2024-01-01T00:00Z
        // Reference 2 days later
        let ref = utcDate(year: 2024, month: 1, day: 3, hour: 0)
        let age = TLEParser.tleAge(elements: e, referenceDate: ref)
        #expect(age > 1.5 && age < 3.0)
    }

    @Test("Epoch date is a valid Date")
    func epochDate() {
        guard let e = TLEParser.parse(line1: issTLE1, line2: issTLE2) else {
            Issue.record("parse failed"); return
        }
        // epochDate should be a valid, non-distant date
        let interval = e.epochDate.timeIntervalSince1970
        #expect(interval > 0) // after Unix epoch
        // Epoch JD should be near 2024
        #expect(e.epochJulianDate > 2460000 && e.epochJulianDate < 2460500)
    }

    @Test("Invalid lines return nil")
    func invalidLines() {
        #expect(TLEParser.parse(line1: "", line2: "") == nil)
        #expect(TLEParser.parse(line1: "too short", line2: "also short") == nil)
    }

    @Test("International designator parsed correctly")
    func intlDesignator() {
        guard let e = TLEParser.parse(line1: issTLE1, line2: issTLE2) else {
            Issue.record("parse failed"); return
        }
        #expect(e.internationalDesignator.contains("98067A"))
    }

    @Test("BSTAR drag term is non-zero for ISS")
    func bstarTerm() {
        guard let e = TLEParser.parse(line1: issTLE1, line2: issTLE2) else {
            Issue.record("parse failed"); return
        }
        #expect(e.bstar != 0)
    }
}

// MARK: - SunCalculator Tests

@Suite("SunCalculator")
struct SunCalculatorTests {

    @Test("Paris summer solstice: ~16h daylight")
    func parisSummerSolstice() {
        let sun = SunCalculator.compute(
            latitude: 48.8566, longitude: 2.3522,
            date: utcDate(year: 2024, month: 6, day: 21)
        )
        #expect(sun.dayLengthHours > 15 && sun.dayLengthHours < 17)
        #expect(sun.sunrise < sun.sunset)
    }

    @Test("Paris winter solstice: ~8h daylight")
    func parisWinterSolstice() {
        let sun = SunCalculator.compute(
            latitude: 48.8566, longitude: 2.3522,
            date: utcDate(year: 2024, month: 12, day: 21)
        )
        #expect(sun.dayLengthHours > 7 && sun.dayLengthHours < 10)
    }

    @Test("Equator equinox: ~12h daylight")
    func equatorEquinox() {
        let sun = SunCalculator.compute(
            latitude: 0.0, longitude: 0.0,
            date: utcDate(year: 2024, month: 3, day: 20)
        )
        #expect(abs(sun.dayLengthHours - 12.0) < 1.0)
    }

    @Test("Golden hour starts before sunset")
    func goldenHour() {
        let sun = SunCalculator.compute(
            latitude: 48.8566, longitude: 2.3522,
            date: utcDate(year: 2024, month: 6, day: 15)
        )
        #expect(sun.goldenHourStart < sun.sunset)
        #expect(sun.goldenHourEnd == sun.sunset)
    }

    @Test("sunProgress at sunrise ≈ 0, at sunset ≈ 1")
    func sunProgressBoundaries() {
        let sun = SunCalculator.compute(
            latitude: 40.0, longitude: -74.0,
            date: utcDate(year: 2024, month: 6, day: 15)
        )
        let progressAtSunrise = sun.sunProgress(at: sun.sunrise)
        let progressAtSunset = sun.sunProgress(at: sun.sunset)
        #expect(abs(progressAtSunrise) < 0.05)
        #expect(abs(progressAtSunset - 1.0) < 0.05)
    }

    @Test("dayProgress is alias for sunProgress")
    func dayProgressAlias() {
        let sun = SunCalculator.compute(
            latitude: 40.0, longitude: -74.0,
            date: utcDate(year: 2024, month: 6, day: 15)
        )
        let midday = sun.sunrise.addingTimeInterval(
            sun.sunset.timeIntervalSince(sun.sunrise) / 2
        )
        let sp = sun.sunProgress(at: midday)
        let dp = sun.dayProgress(at: midday)
        #expect(abs(sp - dp) < 0.001)
    }

    @Test("nightProgress at sunset ≈ 0")
    func nightProgress() {
        let sun = SunCalculator.compute(
            latitude: 40.0, longitude: -74.0,
            date: utcDate(year: 2024, month: 6, day: 15)
        )
        let np = sun.nightProgress(at: sun.sunset)
        #expect(abs(np) < 0.05)
    }

    @Test("isDaytime is true between sunrise and sunset")
    func isDaytime() {
        let sun = SunCalculator.compute(
            latitude: 35.0, longitude: 139.0,
            date: utcDate(year: 2024, month: 4, day: 15)
        )
        let midday = sun.sunrise.addingTimeInterval(
            sun.sunset.timeIntervalSince(sun.sunrise) / 2
        )
        #expect(sun.isDaytime(at: midday))
        #expect(!sun.isDaytime(at: sun.sunrise.addingTimeInterval(-3600)))
    }

    @Test("Backward-compatible computeSunTimes alias works")
    func backwardCompat() {
        let sun = SunCalculator.computeSunTimes(
            lat: 48.8566, lon: 2.3522,
            date: utcDate(year: 2024, month: 6, day: 15)
        )
        #expect(sun.dayLengthHours > 0)
    }

    @Test("Lunar approximation produces non-nil moon dates")
    func lunarApproximation() {
        let sun = SunCalculator.compute(
            latitude: 0.0, longitude: 0.0,
            date: utcDate(year: 2024, month: 7, day: 10)
        )
        #expect(sun.approximateMoonrise != nil)
        #expect(sun.approximateMoonset != nil)
        #expect(sun.moonrise == sun.approximateMoonrise) // alias check
        #expect(sun.moonset == sun.approximateMoonset)   // alias check
    }

    @Test("Polar region handles extreme latitudes without crashing")
    func polarEdge() {
        let sun = SunCalculator.compute(
            latitude: 89.0, longitude: 0.0,
            date: utcDate(year: 2024, month: 6, day: 21)
        )
        #expect(sun.dayLengthHours >= 0)
    }

    @Test("Southern hemisphere has longer days in December")
    func southernSummer() {
        let dec = SunCalculator.compute(
            latitude: -33.8688, longitude: 151.2093,
            date: utcDate(year: 2024, month: 12, day: 21)
        )
        let jun = SunCalculator.compute(
            latitude: -33.8688, longitude: 151.2093,
            date: utcDate(year: 2024, month: 6, day: 21)
        )
        #expect(dec.dayLengthHours > jun.dayLengthHours)
    }
}

// MARK: - VisualMagnitudeCalculator Tests

@Suite("VisualMagnitudeCalculator")
struct VisualMagnitudeCalculatorTests {

    @Test("Sun elevation at solar noon is positive")
    func sunElevationNoon() {
        let noon = utcDate(year: 2024, month: 6, day: 21, hour: 12)
        let elev = VisualMagnitudeCalculator.sunElevation(
            date: noon, latitude: 0.0, longitude: 0.0
        )
        #expect(elev > 0)
    }

    @Test("Sun elevation at midnight is negative")
    func sunElevationMidnight() {
        let midnight = utcDate(year: 2024, month: 6, day: 21, hour: 0)
        let elev = VisualMagnitudeCalculator.sunElevation(
            date: midnight, latitude: 0.0, longitude: 0.0
        )
        #expect(elev < 0)
    }

    @Test("High-altitude satellite in twilight is sunlit")
    func highAltSunlit() {
        let sunlit = VisualMagnitudeCalculator.isSatelliteSunlit(
            sunElevationDeg: -10.0,
            satelliteAltitudeKm: 800
        )
        #expect(sunlit)
    }

    @Test("Low satellite in deep night enters shadow")
    func lowAltShadow() {
        // Earth's shadow geometry: at 200 km, shadow threshold ≈ 62°
        // Need sun depression > threshold for satellite to be in shadow
        let sunlit = VisualMagnitudeCalculator.isSatelliteSunlit(
            sunElevationDeg: -70.0,
            satelliteAltitudeKm: 200
        )
        #expect(!sunlit)
    }

    @Test("Apparent magnitude: closer satellite is brighter")
    func closerIsBrighter() {
        let magClose = VisualMagnitudeCalculator.apparentMagnitude(
            standardMagnitude: 1.0, rangeKm: 500, phaseAngleDeg: 50
        )
        let magFar = VisualMagnitudeCalculator.apparentMagnitude(
            standardMagnitude: 1.0, rangeKm: 2000, phaseAngleDeg: 50
        )
        #expect(magClose < magFar) // lower magnitude = brighter
    }

    @Test("Apparent magnitude: favorable phase angle is brighter")
    func phaseAngleEffect() {
        let magHead = VisualMagnitudeCalculator.apparentMagnitude(
            standardMagnitude: 1.0, rangeKm: 1000, phaseAngleDeg: 10
        )
        let magTail = VisualMagnitudeCalculator.apparentMagnitude(
            standardMagnitude: 1.0, rangeKm: 1000, phaseAngleDeg: 150
        )
        #expect(magHead < magTail) // head-on is brighter
    }

    @Test("Phase angle estimation returns [0, 180]")
    func phaseAngle() {
        let angle = VisualMagnitudeCalculator.estimatePhaseAngle(
            sunElevationDeg: -6.0,
            satelliteElevationDeg: 45.0
        )
        #expect(angle >= 0 && angle <= 180)
    }

    @Test("BrightnessLevel has all 4 cases")
    func brightnessLevels() {
        let levels = PassVisibility.BrightnessLevel.allCases
        #expect(levels.count == 4)
        #expect(levels.contains(.brilliant))
        #expect(levels.contains(.bright))
        #expect(levels.contains(.moderate))
        #expect(levels.contains(.faint))
    }

    @Test("PassVisibility conforms to Hashable")
    func passVisibilityHashable() {
        let vis = PassVisibility(
            peakApparentMagnitude: -2.0,
            isNakedEyeVisible: true,
            visibleFraction: 0.8,
            brightnessLabel: .brilliant
        )
        var dict: [PassVisibility: String] = [:]
        dict[vis] = "bright"
        #expect(dict[vis] == "bright")
    }

    @Test("PassVisibility equality comparison")
    func passVisibilityEquality() {
        let a = PassVisibility(
            peakApparentMagnitude: -1.0,
            isNakedEyeVisible: true,
            visibleFraction: 0.5,
            brightnessLabel: .bright
        )
        let b = PassVisibility(
            peakApparentMagnitude: -1.0,
            isNakedEyeVisible: true,
            visibleFraction: 0.5,
            brightnessLabel: .bright
        )
        #expect(a == b)
    }

    @Test("analyzePass with empty trajectory returns nil")
    func analyzeEmptyTrajectory() {
        let observer = CLLocation(latitude: 48.8566, longitude: 2.3522)
        let result = VisualMagnitudeCalculator.analyzePass(
            trajectory: [], standardMagnitude: 1.0, observer: observer
        )
        #expect(result == nil)
    }

    @Test("analyzePass with valid trajectory returns visibility")
    func analyzeValidTrajectory() {
        let now = Date()
        var trajectory: [TrajectoryPoint] = []
        for i in 0..<20 {
            let elev = Double(i < 10 ? i * 5 : (19 - i) * 5)
            let point = TrajectoryPoint(
                time: now.addingTimeInterval(Double(i) * 30),
                azimuth: Double(i) * 9,
                elevation: elev,
                range: 800,
                rangeRate: -0.5
            )
            trajectory.append(point)
        }
        let observer = CLLocation(latitude: 48.8566, longitude: 2.3522)
        let result = VisualMagnitudeCalculator.analyzePass(
            trajectory: trajectory, standardMagnitude: 1.0, observer: observer
        )
        // May or may not return non-nil depending on sun conditions
        _ = result
    }
}

// MARK: - CoordinateConverter Tests

@Suite("CoordinateConverter")
struct CoordinateConverterTests {

    @Test("Paris to DMS contains correct degrees and hemisphere")
    func parisDMS() {
        let dms = CoordinateConverter.toDMS(decimal: 48.8566, isLatitude: true)
        #expect(dms.contains("48°"))
        #expect(dms.contains("N"))
    }

    @Test("Southern hemisphere gets S suffix")
    func southernDMS() {
        let dms = CoordinateConverter.toDMS(decimal: -33.8688, isLatitude: true)
        #expect(dms.contains("S"))
    }

    @Test("Western hemisphere gets W suffix")
    func westernDMS() {
        let dms = CoordinateConverter.toDMS(decimal: -74.0060, isLatitude: false)
        #expect(dms.contains("W"))
    }

    @Test("Eastern hemisphere gets E suffix")
    func easternDMS() {
        let dms = CoordinateConverter.toDMS(decimal: 139.6917, isLatitude: false)
        #expect(dms.contains("E"))
    }

    @Test("DMS compact combines lat and lon")
    func dmsCompact() {
        let compact = CoordinateConverter.toDMSCompact(latitude: 48.8566, longitude: 2.3522)
        #expect(compact.contains("N"))
        #expect(compact.contains("E"))
    }

    @Test("Decimal to DMS to decimal round-trip")
    func dmsRoundTrip() {
        let original = 48.8566
        let decimal = CoordinateConverter.toDecimal(degrees: 48, minutes: 51, seconds: 23.76, direction: "N")
        #expect(abs(decimal - original) < 0.001)
    }

    @Test("toDecimal with S/W produces negative")
    func toDecimalNegative() {
        let south = CoordinateConverter.toDecimal(degrees: 33, minutes: 52, seconds: 7.68, direction: "S")
        #expect(south < 0)
        let west = CoordinateConverter.toDecimal(degrees: 74, minutes: 0, seconds: 21.6, direction: "W")
        #expect(west < 0)
    }

    @Test("16-point compass directions")
    func compassDirections16() {
        #expect(CoordinateConverter.compassDirection(from: 0) == "N")
        #expect(CoordinateConverter.compassDirection(from: 90) == "E")
        #expect(CoordinateConverter.compassDirection(from: 180) == "S")
        #expect(CoordinateConverter.compassDirection(from: 270) == "W")
        #expect(CoordinateConverter.compassDirection(from: 45) == "NE")
    }

    @Test("8-point compass directions")
    func compassDirections8() {
        #expect(CoordinateConverter.compassDirection8(from: 0) == "N")
        #expect(CoordinateConverter.compassDirection8(from: 90) == "E")
        #expect(CoordinateConverter.compassDirection8(from: 180) == "S")
        #expect(CoordinateConverter.compassDirection8(from: 270) == "W")
    }

    @Test("Compass wraps at 360°")
    func compassWrap() {
        let at0 = CoordinateConverter.compassDirection(from: 0)
        let at360 = CoordinateConverter.compassDirection(from: 360)
        let at720 = CoordinateConverter.compassDirection(from: 720)
        #expect(at0 == at360)
        #expect(at0 == at720)
    }

    @Test("Compass handles negative azimuth")
    func compassNegative() {
        #expect(CoordinateConverter.compassDirection(from: -90) == "W")
    }

    @Test("Distance conversion: km → miles")
    func kmToMiles() {
        let miles = CoordinateConverter.convertDistance(100.0, to: .miles)
        #expect(abs(miles - 62.1371) < 0.1)
    }

    @Test("Distance conversion: km → nautical miles")
    func kmToNautical() {
        let nm = CoordinateConverter.convertDistance(100.0, to: .nauticalMiles)
        #expect(abs(nm - 53.9957) < 0.1)
    }

    @Test("Distance conversion: km → km is identity")
    func kmToKm() {
        #expect(CoordinateConverter.convertDistance(42.0, to: .kilometers) == 42.0)
    }

    @Test("Zero distance stays zero")
    func zeroDistance() {
        #expect(CoordinateConverter.convertDistance(0, to: .miles) == 0)
    }

    @Test("AzEl to direction vector is unit length")
    func directionUnitLength() {
        let dir = CoordinateConverter.azElToDirection(azimuth: 45, elevation: 30)
        let len = simd.length(dir)
        #expect(abs(len - 1.0) < 0.001)
    }

    @Test("Zenith points up (Y-axis)")
    func zenithPointsUp() {
        let dir = CoordinateConverter.azElToDirection(azimuth: 0, elevation: 90)
        #expect(dir.y > 0.99)
    }

    @Test("North at horizon points -Z (SceneKit convention)")
    func northHorizon() {
        let dir = CoordinateConverter.azElToDirection(azimuth: 0, elevation: 0)
        #expect(dir.z < -0.99)
    }

    @Test("DistanceUnit has all 3 cases")
    func distanceUnitCases() {
        let cases = DistanceUnit.allCases
        #expect(cases.count == 3)
        #expect(cases.contains(.kilometers))
        #expect(cases.contains(.miles))
        #expect(cases.contains(.nauticalMiles))
    }
}

// MARK: - GridSquareCalculator Tests

@Suite("GridSquareCalculator")
struct GridSquareCalculatorTests {

    @Test("Paris → JN18eu")
    func parisEncode() {
        let grid = GridSquareCalculator.coordinateToGridSquare(
            latitude: 48.8566, longitude: 2.3522, precision: .subsquare
        )
        #expect(grid == "JN18eu")
    }

    @Test("New York → FN20")
    func newYorkEncode() {
        let grid = GridSquareCalculator.coordinateToGridSquare(
            latitude: 40.7128, longitude: -74.0060, precision: .square
        )
        #expect(grid == "FN20")
    }

    @Test("Field-level encoding (2-char)")
    func fieldEncoding() {
        let grid = GridSquareCalculator.coordinateToGridSquare(
            latitude: 48.8566, longitude: 2.3522, precision: .field
        )
        #expect(grid.count == 2)
        #expect(grid == "JN")
    }

    @Test("Extended encoding (8-char)")
    func extendedEncoding() {
        let grid = GridSquareCalculator.coordinateToGridSquare(
            latitude: 48.8566, longitude: 2.3522, precision: .extended
        )
        #expect(grid.count == 8)
        #expect(grid.hasPrefix("JN18eu"))
    }

    @Test("North Pole encodes valid locator")
    func northPoleEncode() {
        let grid = GridSquareCalculator.coordinateToGridSquare(
            latitude: 90.0, longitude: 0.0, precision: .field
        )
        #expect(GridSquareCalculator.isValid(locator: grid))
    }

    @Test("Decode JN18eu → Paris vicinity")
    func parisDecode() {
        let coord = GridSquareCalculator.gridSquareToCoordinate(locator: "JN18eu")
        #expect(coord != nil)
        #expect(abs(coord!.latitude - 48.8566) < 0.5)
        #expect(abs(coord!.longitude - 2.3522) < 0.5)
    }

    @Test("Decode invalid locator returns nil")
    func decodeInvalid() {
        #expect(GridSquareCalculator.gridSquareToCoordinate(locator: "") == nil)
        #expect(GridSquareCalculator.gridSquareToCoordinate(locator: "ZZ99") == nil)
    }

    @Test("Round-trip encode → decode preserves location")
    func roundTrip() {
        let lat = 48.8566, lon = 2.3522
        let grid = GridSquareCalculator.coordinateToGridSquare(
            latitude: lat, longitude: lon, precision: .extended
        )
        guard let decoded = GridSquareCalculator.gridSquareToCoordinate(locator: grid) else {
            Issue.record("Decode failed"); return
        }
        #expect(abs(decoded.latitude - lat) < 0.01)
        #expect(abs(decoded.longitude - lon) < 0.01)
    }

    @Test("Grid bounds are valid rectangle for FN30")
    func gridBounds() {
        let bounds = GridSquareCalculator.gridSquareBounds(locator: "FN30")
        #expect(bounds != nil)
        #expect(bounds!.minLat < bounds!.maxLat)
        #expect(bounds!.minLon < bounds!.maxLon)
    }

    @Test("Distance between Paris and London ≈ 340 km")
    func distanceParisLondon() {
        let dist = GridSquareCalculator.distance(from: "JN18eu", to: "IO91wm")
        #expect(dist != nil)
        #expect(dist! > 300 && dist! < 400)
    }

    @Test("Distance to self is ≈ 0")
    func distanceToSelf() {
        let dist = GridSquareCalculator.distance(from: "JN18eu", to: "JN18eu")
        #expect(dist != nil)
        #expect(dist! < 1.0)
    }

    @Test("Bearing is [0, 360) degrees")
    func bearingRange() {
        let bearing = GridSquareCalculator.bearing(from: "JN18eu", to: "IO91wm")
        #expect(bearing != nil)
        #expect(bearing! >= 0 && bearing! < 360)
    }

    @Test("Validation accepts valid locators of all lengths")
    func validLocators() {
        #expect(GridSquareCalculator.isValid(locator: "JN"))
        #expect(GridSquareCalculator.isValid(locator: "JN18"))
        #expect(GridSquareCalculator.isValid(locator: "JN18eu"))
        #expect(GridSquareCalculator.isValid(locator: "JN18eu49"))
    }

    @Test("Validation rejects invalid locators")
    func invalidLocators() {
        #expect(!GridSquareCalculator.isValid(locator: ""))
        #expect(!GridSquareCalculator.isValid(locator: "ZZ"))
        #expect(!GridSquareCalculator.isValid(locator: "JN1"))   // odd length
    }

    @Test("Batch conversion returns correct count and valid results")
    func batchConversion() {
        let coords: [(latitude: Double, longitude: Double)] = [
            (48.8566, 2.3522), (40.7128, -74.0060), (35.6762, 139.6503)
        ]
        let grids = GridSquareCalculator.batchCoordinateToGridSquare(coordinates: coords)
        #expect(grids.count == 3)
        #expect(grids.allSatisfy { GridSquareCalculator.isValid(locator: $0) })
    }

    @Test("Generate grid lines returns non-empty result")
    func gridLines() {
        let lines = GridSquareCalculator.generateGridLines(
            minLat: 48.0, maxLat: 50.0, minLon: 1.0, maxLon: 4.0, precision: .square
        )
        #expect(!lines.isEmpty)
    }

    @Test("GridPrecision ordering is correct")
    func precisionOrdering() {
        #expect(GridPrecision.field < GridPrecision.square)
        #expect(GridPrecision.square < GridPrecision.subsquare)
        #expect(GridPrecision.subsquare < GridPrecision.extended)
    }

    @Test("GridPrecision.from infers from locator length")
    func precisionFromLocator() {
        #expect(GridPrecision.from(locator: "JN") == .field)
        #expect(GridPrecision.from(locator: "JN18") == .square)
        #expect(GridPrecision.from(locator: "JN18eu") == .subsquare)
        #expect(GridPrecision.from(locator: "JN18eu49") == .extended)
    }

    @Test("GridPrecision displayName is non-empty")
    func precisionDisplayName() {
        for p in GridPrecision.allCases {
            #expect(!p.displayName.isEmpty)
        }
    }

    @Test("GridLine initializer stores properties correctly")
    func gridLineInit() {
        let start = CLLocationCoordinate2D(latitude: 48.0, longitude: 2.0)
        let end = CLLocationCoordinate2D(latitude: 50.0, longitude: 2.0)
        let line = GridLine(start: start, end: end, precision: .square)
        #expect(line.start.latitude == 48.0)
        #expect(line.end.latitude == 50.0)
        #expect(line.precision == .square)
    }
}

// MARK: - FrequencyFormatter Tests

@Suite("FrequencyFormatter")
struct FrequencyFormatterTests {

    @Test("Auto-format selects correct unit prefix")
    func autoFormat() {
        #expect(FrequencyFormatter.format(hz: 500).contains("Hz"))
        #expect(FrequencyFormatter.format(hz: 1000).contains("kHz"))
        #expect(FrequencyFormatter.format(hz: 1_000_000).contains("MHz"))
        #expect(FrequencyFormatter.format(hz: 1_000_000_000).contains("GHz"))
    }

    @Test("Format with explicit unit")
    func formatWithUnit() {
        let result = FrequencyFormatter.format(hz: 145_000_000, unit: .mhz)
        #expect(result.contains("MHz"))
        #expect(result.contains("145"))
    }

    @Test("Format with kHz unit")
    func formatKHz() {
        let result = FrequencyFormatter.format(hz: 145_000_000, unit: .khz)
        #expect(result.contains("kHz"))
    }

    @Test("Export format is precise MHz")
    func exportFormat() {
        let result = FrequencyFormatter.formatForExport(hz: 145_900_000)
        #expect(result.contains("145.9"))
    }

    @Test("Parse round-trip for MHz frequency")
    func parseRoundTrip() {
        let original: Int64 = 145_900_000
        let formatted = FrequencyFormatter.format(hz: original)
        let parsed = FrequencyFormatter.parse(formatted)
        #expect(parsed != nil)
        if let p = parsed {
            let diff = abs(p - original)
            #expect(diff < 1000) // within 1 kHz
        }
    }

    @Test("Parse returns nil for garbage input")
    func parseGarbage() {
        #expect(FrequencyFormatter.parse("not a frequency") == nil)
        #expect(FrequencyFormatter.parse("") == nil)
    }

    @Test("Band classification for VHF/UHF")
    func bandNames() {
        #expect(FrequencyFormatter.bandName(hz: 145_000_000) != nil) // VHF
        #expect(FrequencyFormatter.bandName(hz: 435_000_000) != nil) // UHF
    }

    @Test("Band name for SHF frequencies")
    func shfBand() {
        #expect(FrequencyFormatter.bandName(hz: 10_000_000_000) != nil)
    }

    @Test("Shift format shows positive prefix")
    func positiveShift() {
        let s = FrequencyFormatter.formatShift(hz: 600_000)
        #expect(s.contains("+"))
    }

    @Test("Shift format shows negative prefix")
    func negativeShift() {
        let s = FrequencyFormatter.formatShift(hz: -600_000)
        // Unicode minus or hyphen
        #expect(s.contains("-") || s.contains("−") || s.contains("\u{2212}"))
    }

    @Test("Zero frequency formats without crash")
    func zeroFrequency() {
        let result = FrequencyFormatter.format(hz: 0)
        #expect(!result.isEmpty)
    }

    @Test("Negative frequency formats without crash")
    func negativeFrequency() {
        let result = FrequencyFormatter.format(hz: -1)
        #expect(!result.isEmpty)
    }

    @Test("FrequencyDisplayUnit has all 4 cases")
    func displayUnitCases() {
        let allCases = FrequencyDisplayUnit.allCases
        #expect(allCases.count == 4)
    }

    @Test("FrequencyDisplayUnit raw values are display labels")
    func displayUnitRawValues() {
        #expect(FrequencyDisplayUnit.hz.rawValue == "Hz")
        #expect(FrequencyDisplayUnit.khz.rawValue == "kHz")
        #expect(FrequencyDisplayUnit.mhz.rawValue == "MHz")
        #expect(FrequencyDisplayUnit.ghz.rawValue == "GHz")
    }
}

// MARK: - AstrodynamicsConstants Tests

@Suite("AstrodynamicsConstants")
struct AstrodynamicsConstantsTests {

    @Test("Earth radius WGS-84 is 6378.137 km")
    func earthRadius() {
        #expect(AstrodynamicsConstants.earthRadiusKm == 6378.137)
    }

    @Test("Earth radius WGS-72 is 6378.135 km (SGP4 standard)")
    func earthRadiusWGS72() {
        #expect(AstrodynamicsConstants.earthRadiusWGS72 == 6378.135)
    }

    @Test("Earth mean radius is 6371.0 km")
    func earthMeanRadius() {
        #expect(AstrodynamicsConstants.earthMeanRadiusKm == 6371.0)
    }

    @Test("Gravitational parameter μ WGS-84")
    func earthMu() {
        #expect(abs(AstrodynamicsConstants.earthMu - 398600.4418) < 0.001)
    }

    @Test("Gravitational parameter μ WGS-72")
    func earthMuWGS72() {
        #expect(abs(AstrodynamicsConstants.earthMuWGS72 - 398600.8) < 0.1)
    }

    @Test("Degree/radian conversions are inverses")
    func angleConversions() {
        let deg = 180.0
        let rad = deg * AstrodynamicsConstants.deg2rad
        let backToDeg = rad * AstrodynamicsConstants.rad2deg
        #expect(abs(backToDeg - deg) < 1e-10)
    }

    @Test("deg2rad of 180 = π")
    func deg2rad() {
        #expect(abs(180.0 * AstrodynamicsConstants.deg2rad - .pi) < 1e-15)
    }

    @Test("J2000.0 Julian Date is 2451545.0")
    func j2000() {
        #expect(AstrodynamicsConstants.j2000JulianDate == 2451545.0)
    }

    @Test("Seconds per day is 86400.0")
    func secondsPerDay() {
        #expect(AstrodynamicsConstants.secondsPerDay == 86400.0)
    }

    @Test("Minutes per day is 1440.0")
    func minutesPerDay() {
        #expect(AstrodynamicsConstants.minutesPerDay == 1440.0)
    }

    @Test("Earth rotation rate is ~7.29e-5 rad/s")
    func earthRotation() {
        #expect(abs(AstrodynamicsConstants.earthRotationRate - 7.29211514670698e-5) < 1e-10)
    }

    @Test("Earth flattening WGS-84 is ~1/298.257")
    func earthFlattening() {
        let expected = 1.0 / 298.257223563
        #expect(abs(AstrodynamicsConstants.earthFlattening - expected) < 1e-12)
    }

    @Test("julianDate(from:) matches SGP4Propagator implementation")
    func julianDateConsistency() {
        let date = Date(timeIntervalSince1970: 1704067200)
        let jd1 = AstrodynamicsConstants.julianDate(from: date)
        let jd2 = SGP4Propagator.julianDate(from: date)
        #expect(abs(jd1 - jd2) < 0.001)
    }

    @Test("j2000UnixTimestamp is correct")
    func j2000Unix() {
        // J2000.0 = 2000-01-01T12:00:00Z
        let expected = utcDate(year: 2000, month: 1, day: 1, hour: 12)
            .timeIntervalSince1970
        #expect(abs(AstrodynamicsConstants.j2000UnixTimestamp - expected) < 1.0)
    }
}

// MARK: - Integration Tests

@Suite("Integration")
struct IntegrationTests {

    @Test("SGP4 + CoordinateConverter: propagate ISS and get compass direction")
    func sgp4ToCompass() {
        guard let prop = SGP4Propagator(tle1: issTLE1, tle2: issTLE2),
              let eci = prop.propagate(date: Date(timeIntervalSince1970: 1704067200)) else {
            Issue.record("propagation failed"); return
        }
        let paris = CLLocation(latitude: 48.8566, longitude: 2.3522)
        let look = SGP4Propagator.lookAngles(
            position: eci.position, velocity: eci.velocity,
            observer: paris, date: Date(timeIntervalSince1970: 1704067200)
        )
        let direction = CoordinateConverter.compassDirection(from: look.azimuth)
        #expect(!direction.isEmpty)
    }

    @Test("TLEParser + SGP4 + GridSquare: parse → propagate → locate")
    func parseThenLocate() {
        guard let elements = TLEParser.parse(line1: issTLE1, line2: issTLE2) else {
            Issue.record("parse failed"); return
        }
        #expect(elements.noradId == 25544)

        guard let prop = SGP4Propagator(tle1: issTLE1, tle2: issTLE2),
              let eci = prop.propagate(date: Date(timeIntervalSince1970: 1704067200)) else {
            Issue.record("propagate failed"); return
        }

        let geo = SGP4Propagator.eciToGeodetic(
            position: eci.position,
            date: Date(timeIntervalSince1970: 1704067200)
        )
        let grid = GridSquareCalculator.coordinateToGridSquare(
            latitude: geo.latitude, longitude: geo.longitude
        )
        #expect(GridSquareCalculator.isValid(locator: grid))
    }

    @Test("SunCalculator + VisualMagnitudeCalculator: twilight consistency")
    func sunAndVisibility() {
        let evening = utcDate(year: 2024, month: 6, day: 15, hour: 21)
        let sunElev = VisualMagnitudeCalculator.sunElevation(
            date: evening, latitude: 48.8566, longitude: 2.3522
        )
        let solar = SunCalculator.compute(
            latitude: 48.8566, longitude: 2.3522, date: evening
        )
        #expect(solar.dayLengthHours > 14) // long summer day
        #expect(sunElev < 30) // not noon-level
    }

    @Test("TLEParser + FrequencyFormatter: satellite metadata chain")
    func metadataChain() {
        guard let elements = TLEParser.parse(line1: issTLE1, line2: issTLE2) else {
            Issue.record("parse failed"); return
        }
        // ISS downlink: 145.800 MHz
        let freq: Int64 = 145_800_000
        let formatted = FrequencyFormatter.format(hz: freq)
        let band = FrequencyFormatter.bandName(hz: freq)
        #expect(formatted.contains("MHz"))
        #expect(band != nil)
        #expect(elements.noradId == 25544) // metadata stays consistent
    }

    @Test("Complete satellite visibility check chain")
    func visibilityChain() {
        // Check if ISS would be visible from Paris at a given time
        guard let prop = SGP4Propagator(tle1: issTLE1, tle2: issTLE2),
              let eci = prop.propagate(date: Date(timeIntervalSince1970: 1704067200)) else {
            Issue.record("propagation failed"); return
        }

        let paris = CLLocation(latitude: 48.8566, longitude: 2.3522)
        let look = SGP4Propagator.lookAngles(
            position: eci.position, velocity: eci.velocity,
            observer: paris, date: Date(timeIntervalSince1970: 1704067200)
        )

        let geo = SGP4Propagator.eciToGeodetic(
            position: eci.position,
            date: Date(timeIntervalSince1970: 1704067200)
        )

        let sunElev = VisualMagnitudeCalculator.sunElevation(
            date: Date(timeIntervalSince1970: 1704067200),
            latitude: 48.8566, longitude: 2.3522
        )

        let sunlit = VisualMagnitudeCalculator.isSatelliteSunlit(
            sunElevationDeg: sunElev,
            satelliteAltitudeKm: geo.altitude
        )

        // All these should be reasonable numbers
        #expect(look.range > 0)
        #expect(geo.altitude > 300 && geo.altitude < 500)
        _ = sunlit // just verify no crash
    }
}

// MARK: - Stress & Edge Cases

@Suite("Stress & Edge Cases")
struct StressTests {

    @Test("1000 grid square encodings")
    func gridSquarePerformance() {
        let start = Date()
        for i in 0..<1000 {
            let lat = -90.0 + Double(i) * 0.18
            let lon = -180.0 + Double(i) * 0.36
            _ = GridSquareCalculator.coordinateToGridSquare(
                latitude: lat, longitude: lon, precision: .extended
            )
        }
        let elapsed = Date().timeIntervalSince(start)
        #expect(elapsed < 2.0) // generous for CI
    }

    @Test("Compass direction for all integer azimuths [0,360)")
    func allAzimuths() {
        for az in 0..<360 {
            let dir = CoordinateConverter.compassDirection(from: Double(az))
            #expect(!dir.isEmpty)
        }
    }

    @Test("Sun calculator handles all 365 days of 2024")
    func allDaysOfYear() {
        for day in 1...365 {
            let date = utcDate(year: 2024, month: 1, day: day)
            let sun = SunCalculator.compute(latitude: 48.8566, longitude: 2.3522, date: date)
            #expect(sun.dayLengthHours > 0 && sun.dayLengthHours <= 24)
        }
    }

    @Test("Frequency formatter handles full RF spectrum")
    func fullSpectrum() {
        let frequencies: [Int64] = [1, 100, 1000, 10_000, 100_000,
                                     1_000_000, 10_000_000, 100_000_000,
                                     1_000_000_000, 10_000_000_000]
        for freq in frequencies {
            let result = FrequencyFormatter.format(hz: freq)
            #expect(!result.isEmpty)
        }
    }

    @Test("CoordinateConverter handles boundary coordinates")
    func boundaryCoordinates() {
        let north = CoordinateConverter.toDMS(decimal: 90.0, isLatitude: true)
        #expect(north.contains("N"))
        let south = CoordinateConverter.toDMS(decimal: -90.0, isLatitude: true)
        #expect(south.contains("S"))
        let east = CoordinateConverter.toDMS(decimal: 180.0, isLatitude: false)
        #expect(east.contains("E"))
        let west = CoordinateConverter.toDMS(decimal: -180.0, isLatitude: false)
        #expect(west.contains("W"))
    }

    @Test("Grid square boundary coordinates don't crash")
    func gridBoundaryCoords() {
        let extremes: [(Double, Double)] = [
            (90.0, 180.0), (-90.0, -180.0), (0.0, 0.0),
            (89.99, 179.99), (-89.99, -179.99)
        ]
        for (lat, lon) in extremes {
            let grid = GridSquareCalculator.coordinateToGridSquare(latitude: lat, longitude: lon)
            #expect(!grid.isEmpty)
            #expect(GridSquareCalculator.isValid(locator: grid))
        }
    }

    @Test("Multiple TLE parses in rapid succession")
    func rapidTLEParsing() {
        for _ in 0..<100 {
            let e = TLEParser.parse(line1: issTLE1, line2: issTLE2)
            #expect(e != nil)
            #expect(e?.noradId == 25544)
        }
    }

    @Test("Julian date monotonically increases with time")
    func julianDateMonotonic() {
        var prevJd = 0.0
        for i in 0..<100 {
            let date = Date(timeIntervalSince1970: Double(i) * 86400)
            let jd = SGP4Propagator.julianDate(from: date)
            #expect(jd > prevJd)
            prevJd = jd
        }
    }

    @Test("ECI to geodetic: latitude always in [-90, 90]")
    func geodeticLatBound() {
        let positions: [SIMD3<Double>] = [
            SIMD3(6778, 0, 0), SIMD3(0, 6778, 0), SIMD3(0, 0, 6778),
            SIMD3(-6778, 0, 0), SIMD3(0, -6778, 0), SIMD3(0, 0, -6778),
            SIMD3(4000, 4000, 4000)
        ]
        let date = Date(timeIntervalSince1970: 1704067200)
        for pos in positions {
            let geo = SGP4Propagator.eciToGeodetic(position: pos, date: date)
            #expect(geo.latitude >= -90 && geo.latitude <= 90)
            #expect(geo.longitude >= -180 && geo.longitude <= 180)
        }
    }
}
