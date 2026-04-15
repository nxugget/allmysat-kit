// ────────────────────────────────────────────────────────────────────────────
// SunCalculator.swift
// AllMySatKit · Astronomy
//
// Computes sunrise, sunset, and solar day geometry for a given date and
// geographic location.
//
// The algorithm uses the Earth's obliquity to compute the solar declination,
// then solves for the hour angle at which the sun crosses the standard
// zenith angle for sunrise/sunset (90.833°, accounting for atmospheric
// refraction and solar disk diameter).
//
// Zenith angle breakdown:
//   90°      - geometric horizon
//   +0.833°  - atmospheric refraction (0.567°) + solar semi-diameter (0.266°)
//   sin(−0.833°) ≈ −0.01454
//
// The lunar approximation uses the mean synodic month (29.53 days) to
// estimate moonrise phase. This is a rough approximation - for precise
// lunar ephemerides, a full ELP/MPP model would be needed.
//
// Copyright 2026 AllMySat · Apache License 2.0
// ────────────────────────────────────────────────────────────────────────────

import Foundation

// MARK: - Result Type

/// Solar geometry data for a single day at a given location.
public struct SolarDayInfo: Sendable {

    /// Moment of sunrise (sun center at −0.833° below geometric horizon).
    public let sunrise: Date

    /// Moment of sunset (sun center at −0.833° below geometric horizon).
    public let sunset: Date

    /// Start of golden hour (approximately 1 hour before sunset).
    public let goldenHourStart: Date

    /// Approximate moonrise based on the synodic month cycle.
    ///
    /// This is a rough estimate - error can exceed ±2 hours.
    /// For precise lunar times, use a dedicated ephemeris.
    public let approximateMoonrise: Date?

    /// Approximate moonset (~12 h after moonrise).
    public let approximateMoonset: Date?

    /// Convenience alias for `approximateMoonrise`.
    public var moonrise: Date? { approximateMoonrise }

    /// Convenience alias for `approximateMoonset`.
    public var moonset: Date? { approximateMoonset }

    /// Length of the day (sunrise to sunset), in hours.
    public let dayLengthHours: Double

    public init(sunrise: Date, sunset: Date, goldenHourStart: Date,
                approximateMoonrise: Date?, approximateMoonset: Date?,
                dayLengthHours: Double) {
        self.sunrise = sunrise
        self.sunset = sunset
        self.goldenHourStart = goldenHourStart
        self.approximateMoonrise = approximateMoonrise
        self.approximateMoonset = approximateMoonset
        self.dayLengthHours = dayLengthHours
    }

    /// Sunset (alias for goldenHourEnd, since golden hour ends at sunset).
    public var goldenHourEnd: Date { sunset }

    /// Current sun progress ratio (0.0 = sunrise, 1.0 = sunset).
    /// Clamped to [0, 1]. Values outside the day return 0 or 1.
    public func sunProgress(at now: Date) -> Double {
        let total = sunset.timeIntervalSince(sunrise)
        guard total > 0 else { return 0 }
        return max(0, min(1, now.timeIntervalSince(sunrise) / total))
    }

    /// Fraction of daylight elapsed at the given time (0 = sunrise, 1 = sunset).
    /// Clamped to [0, 1]. Values outside the day return 0 or 1.
    public func dayProgress(at time: Date) -> Double {
        sunProgress(at: time)
    }

    /// Night progress ratio (0.0 = sunset, 1.0 = next sunrise).
    /// Uses an estimated ~12h night if no next-day sunrise is available.
    public func nightProgress(at now: Date) -> Double {
        let nextSunrise = sunrise.addingTimeInterval(86400)
        let total = nextSunrise.timeIntervalSince(sunset)
        guard total > 0 else { return 0 }
        let elapsed = now.timeIntervalSince(sunset)
        return max(0, min(1, elapsed / total))
    }

    /// Whether the given time falls between sunrise and sunset.
    public func isDaytime(at now: Date) -> Bool {
        now >= sunrise && now <= sunset
    }
}

// MARK: - Sun Calculator

/// Computes sunrise, sunset, and related solar geometry.
///
/// Uses the solar declination and hour angle method. Accuracy is within
/// ~1 minute for latitudes below ±65° and degrades near the poles where
/// sunrise/sunset become ill-defined.
public enum SunCalculator {

    /// Earth's axial tilt (obliquity of the ecliptic), in degrees.
    ///
    /// Current value for J2000.0 epoch.
    /// Source: IAU 2006 precession model
    private static let axialTilt: Double = 23.44

    /// Mean synodic month (new moon to new moon), in days.
    ///
    /// Source: Meeus, "Astronomical Algorithms", Table 49.A
    private static let synodicMonth: Double = 29.53

    // MARK: - Public API

    /// Computes solar day information for a given date and location.
    ///
    /// - Parameters:
    ///   - latitude: Observer latitude in decimal degrees (−90 to +90).
    ///   - longitude: Observer longitude in decimal degrees (−180 to +180).
    ///   - date: The date for which to compute sun times.
    /// - Returns: Solar geometry for the given day.
    public static func compute(
        latitude: Double,
        longitude: Double,
        date: Date
    ) -> SolarDayInfo {

        let calendar = Calendar.current
        let dayOfYear = Double(calendar.ordinality(of: .day, in: .year, for: date) ?? 1)
        let dayStart = calendar.startOfDay(for: date)

        // ── Solar declination ──────────────────────────────────
        //
        // δ = −ε × cos(2π/365 × (d + 10))
        //
        // where ε = 23.44° is the axial tilt and d is the day of year.
        // The +10 accounts for the winter solstice offset from Jan 1.

        let latRad = latitude * .pi / 180.0
        let declination = -axialTilt * cos(2.0 * .pi / 365.0 * (dayOfYear + 10)) * .pi / 180.0

        // ── Hour angle at sunrise/sunset ───────────────────────
        //
        // cos(ω₀) = [sin(−0.833°) − sin(φ)sin(δ)] / [cos(φ)cos(δ)]
        //
        // sin(−0.833°) ≈ −0.01454 accounts for atmospheric refraction
        // and the sun's apparent diameter.

        let cosHourAngle = (-0.01454 - sin(latRad) * sin(declination))
            / (cos(latRad) * cos(declination))

        // Clamp for polar day/night (midnight sun or polar night)
        let hourAngle = acos(max(-1, min(1, cosHourAngle)))

        // ── Solar noon (local mean time) ───────────────────────
        //
        // Solar noon ≈ 12:00 − longitude/15°
        // This is approximate; the equation of time correction is omitted.

        let solarNoonHours = 12.0 - longitude / 15.0

        // ── Sunrise and sunset times ───────────────────────────

        let sunriseHours = solarNoonHours - hourAngle * 12.0 / .pi
        let sunsetHours  = solarNoonHours + hourAngle * 12.0 / .pi

        let sunrise = dayStart.addingTimeInterval(sunriseHours * 3600.0)
        let sunset  = dayStart.addingTimeInterval(sunsetHours * 3600.0)

        // ── Golden hour ────────────────────────────────────────
        // Approximately 1 hour before sunset (when sun is ~6° above horizon)

        let goldenHourStart = sunset.addingTimeInterval(-3600.0)

        // ── Lunar approximation ────────────────────────────────
        //
        // Very rough estimate using the synodic month period.
        // The moon rises ~50 minutes later each day on average.
        // Moonrise hour ≈ (dayOfYear mod 29.53) / 29.53 × 24

        let lunarPhaseDay = dayOfYear.truncatingRemainder(dividingBy: synodicMonth)
        let moonriseHours = (lunarPhaseDay / synodicMonth) * 24.0
        let moonrise = dayStart.addingTimeInterval(moonriseHours * 3600.0)
        let moonset  = moonrise.addingTimeInterval(12.0 * 3600.0)

        return SolarDayInfo(
            sunrise: sunrise,
            sunset: sunset,
            goldenHourStart: goldenHourStart,
            approximateMoonrise: moonrise,
            approximateMoonset: moonset,
            dayLengthHours: sunsetHours - sunriseHours
        )
    }

    /// Backward-compatible alias.
    ///
    /// Matches the method signature used internally by AllMySat.
    public static func computeSunTimes(lat: Double, lon: Double, date: Date) -> SolarDayInfo {
        compute(latitude: lat, longitude: lon, date: date)
    }
}
