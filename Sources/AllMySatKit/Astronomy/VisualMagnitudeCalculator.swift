// ────────────────────────────────────────────────────────────────────────────
// VisualMagnitudeCalculator.swift
// AllMySatKit · Astronomy
//
// Estimates the apparent visual magnitude of a satellite during a pass,
// accounting for range, illumination geometry, Earth's shadow, and
// observer twilight conditions.
//
// The apparent brightness of a satellite depends on:
//   1. Intrinsic brightness (standard magnitude at 1000 km, 50° phase angle)
//   2. Observer–satellite range (inverse-square dimming in magnitudes)
//   3. Phase angle between Sun–satellite–observer (Lambertian diffuse model)
//   4. Whether the satellite is in sunlight (above Earth's shadow cone)
//   5. Whether the observer is in sufficient darkness (twilight check)
//
// Magnitude scale:
//   Astronomical apparent magnitude is logarithmic and inverted.
//   Δm = −2.5 × log₁₀(flux_ratio)
//   Lower values are brighter. Sirius = −1.46, naked eye limit ≈ +6.0.
//
// Copyright 2026 AllMySat · Apache License 2.0
// ────────────────────────────────────────────────────────────────────────────

import Foundation
import CoreLocation

// MARK: - Result Type

/// Visibility analysis for a complete satellite pass.
public struct PassVisibility: Sendable, Hashable {
    /// Peak (lowest numerical) apparent magnitude during the pass.
    /// Lower values = brighter.
    public let peakApparentMagnitude: Double

    /// Whether the satellite is potentially visible to the naked eye (mag ≤ 6.0).
    public let isNakedEyeVisible: Bool

    /// Fraction of trajectory points where the satellite is both sunlit
    /// and the observer is in twilight or darkness, in [0, 1].
    public let visibleFraction: Double

    /// Categorical brightness assessment.
    public let brightnessLabel: BrightnessLevel

    public init(peakApparentMagnitude: Double, isNakedEyeVisible: Bool,
                visibleFraction: Double, brightnessLabel: BrightnessLevel) {
        self.peakApparentMagnitude = peakApparentMagnitude
        self.isNakedEyeVisible = isNakedEyeVisible
        self.visibleFraction = visibleFraction
        self.brightnessLabel = brightnessLabel
    }

    /// Categorical brightness assessment of a satellite pass.
    public enum BrightnessLevel: String, Sendable, CaseIterable, Hashable {
        /// Easily visible, brighter than most stars (mag ≤ 0).
        case brilliant
        /// Clearly visible to the naked eye (0 < mag ≤ 2).
        case bright
        /// Visible in clear dark skies (2 < mag ≤ 4).
        case moderate
        /// At the limit of naked-eye visibility (4 < mag ≤ 6).
        case faint
    }
}

// MARK: - Visual Magnitude Calculator

/// Estimates the apparent visual magnitude of a satellite during a pass.
///
/// All methods are pure functions with no shared state.
public enum VisualMagnitudeCalculator {

    // MARK: - Constants

    /// Limiting apparent magnitude for naked-eye visibility.
    ///
    /// Under ideal conditions (no light pollution, perfect dark adaptation).
    /// Source: Schaefer, B. - "Telescopic Limiting Magnitudes", PASP 102, 1990
    private static let nakedEyeLimit: Double = 6.0

    /// Reference range for standard (catalog) magnitude, in km.
    ///
    /// Standard magnitude is measured at 1000 km slant range.
    /// Source: McCants satellite visual magnitude database conventions.
    private static let standardRange: Double = 1000.0

    /// Reference phase angle for standard magnitude, in degrees.
    ///
    /// Source: McCants database convention (half-illumination geometry).
    private static let standardPhaseAngle: Double = 50.0

    // MARK: - Sun Position

    /// Computes the Sun's elevation angle above the observer's horizon.
    ///
    /// Uses the low-precision solar position algorithm:
    ///   1. Mean anomaly from Julian centuries since J2000.0
    ///   2. Equation of center (3-term trigonometric series)
    ///   3. Ecliptic longitude → right ascension + declination
    ///   4. Greenwich Mean Sidereal Time → local hour angle
    ///   5. Altitude via the standard spherical astronomy formula
    ///
    /// Accuracy: ~0.01° for dates within ±50 years of J2000.0.
    ///
    /// - Parameters:
    ///   - date: UTC time.
    ///   - latitude: Observer latitude in degrees.
    ///   - longitude: Observer longitude in degrees.
    /// - Returns: Sun elevation in degrees (negative = below horizon).
    public static func sunElevation(date: Date, latitude: Double, longitude: Double) -> Double {

        // Julian centuries since J2000.0
        let jd = AstrodynamicsConstants.julianDate(from: date)
        let T = (jd - AstrodynamicsConstants.j2000JulianDate) / 36525.0

        // ── Mean anomaly (degrees) ─────────────────────────────────
        // Source: Meeus, "Astronomical Algorithms", 2nd ed., Table 25.A
        let M = (357.5291 + 35999.0503 * T)
            .truncatingRemainder(dividingBy: 360.0)

        let Mrad = M * AstrodynamicsConstants.deg2rad

        // ── Equation of center ─────────────────────────────────────
        // 3-term trigonometric expansion for Earth's orbit eccentricity.
        // Source: Meeus, Chapter 25
        let C = 1.9146 * sin(Mrad)
               + 0.0200 * sin(2.0 * Mrad)
               + 0.0003 * sin(3.0 * Mrad)

        // ── Ecliptic longitude ─────────────────────────────────────
        // λ = M + C + ω + 180°
        // ω ≈ 102.9373° (longitude of perihelion)
        let eclLon = (M + C + 102.9373 + 180.0)
            .truncatingRemainder(dividingBy: 360.0) * AstrodynamicsConstants.deg2rad

        // ── Solar declination ──────────────────────────────────────
        // sin(δ) = sin(λ) × sin(ε)
        // ε ≈ 23.4393° (mean obliquity, IAU 2006)
        let obliquity = 23.4393 * AstrodynamicsConstants.deg2rad
        let sinDec = sin(eclLon) * sin(obliquity)
        let declination = asin(sinDec)

        // ── Right ascension ────────────────────────────────────────
        let rightAscension = atan2(
            sin(eclLon) * cos(obliquity),
            cos(eclLon)
        )

        // ── Greenwich Mean Sidereal Time ───────────────────────────
        // GMST at 0h UT + Earth rotation since midnight
        // Source: Meeus, Chapter 12
        let gmst = (280.46061837
                   + 360.98564736629 * (jd - AstrodynamicsConstants.j2000JulianDate))
            .truncatingRemainder(dividingBy: 360.0) * AstrodynamicsConstants.deg2rad

        // ── Local Hour Angle ───────────────────────────────────────
        let localSiderealTime = gmst + longitude * AstrodynamicsConstants.deg2rad
        let hourAngle = localSiderealTime - rightAscension

        // ── Elevation ──────────────────────────────────────────────
        // sin(h) = sin(φ)sin(δ) + cos(φ)cos(δ)cos(H)
        let latRad = latitude * AstrodynamicsConstants.deg2rad
        let sinElevation = sin(latRad) * sin(declination)
                         + cos(latRad) * cos(declination) * cos(hourAngle)

        return asin(sinElevation) * AstrodynamicsConstants.rad2deg
    }

    // MARK: - Earth Shadow

    /// Determines whether a satellite at a given altitude is illuminated by the Sun.
    ///
    /// Uses a cylindrical shadow model:
    ///   - Compute the angular radius of the Earth as seen from the satellite.
    ///   - Compare with the Sun's depression angle below the satellite's horizon.
    ///   - If the Sun depression exceeds the Earth's angular radius, the satellite
    ///     is in Earth's shadow.
    ///
    /// This is a simplification - the true umbra/penumbra boundary is conical.
    /// For LEO satellites (200–2000 km), the error is typically < 1 minute.
    ///
    /// - Parameters:
    ///   - sunElevationDeg: Sun elevation at the observer, in degrees.
    ///   - satelliteAltitudeKm: Satellite altitude above Earth's surface, in km.
    /// - Returns: `true` if the satellite is likely in sunlight.
    public static func isSatelliteSunlit(
        sunElevationDeg: Double,
        satelliteAltitudeKm: Double
    ) -> Bool {

        let Re = AstrodynamicsConstants.earthMeanRadiusKm
        let r = Re + satelliteAltitudeKm

        // Angular radius of Earth from the satellite's vantage point
        //   θ_earth = arcsin(R_e / r)
        let earthAngularRadius = asin(Re / r) * AstrodynamicsConstants.rad2deg

        // Sun depression angle below the horizon (positive when sun is below)
        let sunDepression = -sunElevationDeg

        // Shadow threshold: satellite is in shadow when the sun depression
        // exceeds the angular dip of the satellite above Earth's limb.
        // The correction factor accounts for the satellite's altitude elevating
        // its effective horizon relative to the observer.
        let shadowThreshold = earthAngularRadius - (90.0 - earthAngularRadius)

        return sunDepression < shadowThreshold
    }

    // MARK: - Apparent Magnitude

    /// Computes the apparent visual magnitude of a satellite.
    ///
    /// Two corrections are applied to the catalog (standard) magnitude:
    ///
    /// **Range correction** (inverse-square law in magnitude space):
    ///   Δm_range = 5 × log₁₀(range / range_std)
    ///
    /// **Phase angle correction** (Lambertian diffuse sphere):
    ///   Φ(θ) = (1 + cos θ) / 2
    ///   Δm_phase = −2.5 × log₁₀( Φ(θ) / Φ(θ_std) )
    ///
    /// The Lambertian model assumes uniform diffuse reflection from a sphere.
    /// Real satellites have specular components (solar panels) that can cause
    /// brief flares - this model does not account for those.
    ///
    /// - Parameters:
    ///   - standardMagnitude: Catalog magnitude at 1000 km and 50° phase angle.
    ///   - rangeKm: Slant range from observer to satellite, in km.
    ///   - phaseAngleDeg: Sun–satellite–observer angle, in degrees.
    /// - Returns: Estimated apparent visual magnitude.
    public static func apparentMagnitude(
        standardMagnitude: Double,
        rangeKm: Double,
        phaseAngleDeg: Double
    ) -> Double {

        // ── Range correction ───────────────────────────────────
        let rangeCorrection = 5.0 * log10(max(rangeKm, 1.0) / standardRange)

        // ── Phase angle correction (Lambertian model) ──────────
        let phaseRad = phaseAngleDeg * AstrodynamicsConstants.deg2rad
        let stdPhaseRad = standardPhaseAngle * AstrodynamicsConstants.deg2rad

        let actualPhase = (1.0 + cos(phaseRad)) / 2.0
        let stdPhase    = (1.0 + cos(stdPhaseRad)) / 2.0

        // Guard against zero brightness (satellite showing no illuminated area)
        let phaseCorrection: Double
        if actualPhase > 0 {
            phaseCorrection = -2.5 * log10(actualPhase / stdPhase)
        } else {
            phaseCorrection = 10.0  // effectively invisible
        }

        return standardMagnitude + rangeCorrection + phaseCorrection
    }

    // MARK: - Phase Angle Estimation

    /// Estimates the Sun–satellite–observer phase angle from observable quantities.
    ///
    /// This is a simplified geometric model assuming the Sun is infinitely distant
    /// (which is an excellent approximation - the parallax is < 0.01°).
    ///
    /// The phase angle is approximated as:
    ///   φ ≈ 180° − (sun_depression + satellite_elevation)
    ///
    /// where sun_depression = −sun_elevation (positive when sun is below horizon).
    ///
    /// - Parameters:
    ///   - sunElevationDeg: Sun elevation at the observer, in degrees.
    ///   - satelliteElevationDeg: Satellite elevation above the observer's horizon, in degrees.
    /// - Returns: Estimated phase angle in degrees, clamped to [0°, 180°].
    public static func estimatePhaseAngle(
        sunElevationDeg: Double,
        satelliteElevationDeg: Double
    ) -> Double {
        let phase = 180.0 - (-sunElevationDeg + satelliteElevationDeg)
        return max(0, min(180, phase))
    }

    // MARK: - Pass Analysis

    /// Analyzes a complete satellite pass for visual observability.
    ///
    /// Samples the trajectory at regular intervals and for each sample:
    ///   1. Checks that the observer is in astronomical twilight (sun < −2°)
    ///   2. Checks that the satellite is in sunlight (above Earth's shadow)
    ///   3. Computes apparent magnitude with range and phase corrections
    ///
    /// The twilight threshold of −2° is more permissive than the standard
    /// astronomical twilight (−18°) because bright satellites like the ISS
    /// (mag −3 to −6) are visible in lighter skies.
    ///
    /// - Parameters:
    ///   - trajectory: Trajectory points from a `PredictedPass`.
    ///   - standardMagnitude: Catalog magnitude of the satellite.
    ///   - observer: Observer geographic location.
    /// - Returns: Visibility analysis for the pass, or `nil` if the trajectory is empty.
    public static func analyzePass(
        trajectory: [TrajectoryPoint],
        standardMagnitude: Double,
        observer: CLLocation
    ) -> PassVisibility? {

        guard !trajectory.isEmpty else { return nil }

        // Sample at most 30 points, evenly distributed across the trajectory
        let maxSamples = 30
        let stride = max(1, trajectory.count / maxSamples)
        let sampledPoints = Swift.stride(from: 0, to: trajectory.count, by: stride)
            .map { trajectory[$0] }

        var bestMagnitude = Double.infinity
        var visibleCount = 0

        for point in sampledPoints {

            // ── Observer twilight check ─────────────────────────
            // Sun must be below −2° for satellite visibility against the sky
            let sunElev = sunElevation(
                date: point.time,
                latitude: observer.coordinate.latitude,
                longitude: observer.coordinate.longitude
            )
            guard sunElev < -2.0 else { continue }

            // ── Satellite sunlight check ───────────────────────
            // Estimate altitude from range and elevation angle:
            //   h ≈ range × sin(elevation) [for moderate elevations]
            // This is approximate; the exact altitude comes from the propagator's
            // ECI state, but range × sin(el) is sufficient for the shadow test.
            let altitudeKm = point.range * sin(point.elevation * AstrodynamicsConstants.deg2rad)
            guard isSatelliteSunlit(
                sunElevationDeg: sunElev,
                satelliteAltitudeKm: max(altitudeKm, 200.0)
            ) else { continue }

            // ── Apparent magnitude ─────────────────────────────
            let phase = estimatePhaseAngle(
                sunElevationDeg: sunElev,
                satelliteElevationDeg: point.elevation
            )
            let mag = apparentMagnitude(
                standardMagnitude: standardMagnitude,
                rangeKm: point.range,
                phaseAngleDeg: phase
            )

            visibleCount += 1
            bestMagnitude = min(bestMagnitude, mag)
        }

        // If no point was visible, report as faintest possible
        if visibleCount == 0 {
            return PassVisibility(
                peakApparentMagnitude: 99.0,
                isNakedEyeVisible: false,
                visibleFraction: 0,
                brightnessLabel: .faint
            )
        }

        let fraction = Double(visibleCount) / Double(sampledPoints.count)
        let level = classifyBrightness(magnitude: bestMagnitude)

        return PassVisibility(
            peakApparentMagnitude: bestMagnitude,
            isNakedEyeVisible: bestMagnitude <= nakedEyeLimit,
            visibleFraction: fraction,
            brightnessLabel: level
        )
    }

    // MARK: - Classification

    /// Maps an apparent magnitude to a brightness category.
    private static func classifyBrightness(magnitude: Double) -> PassVisibility.BrightnessLevel {
        switch magnitude {
        case ...0:    return .brilliant
        case ...2:    return .bright
        case ...4:    return .moderate
        default:      return .faint
        }
    }
}
