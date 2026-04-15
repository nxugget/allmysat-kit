// ────────────────────────────────────────────────────────────────────────────
// AstrodynamicsConstants.swift
// AllMySatKit
//
// Fundamental physical and astronomical constants used across the library.
// Each constant is annotated with its source and standard.
//
// Standards used:
//   WGS-72  - World Geodetic System 1972 (SGP4 standard per Spacetrack Report #3)
//   WGS-84  - World Geodetic System 1984 (modern geodetic reference)
//   IAU     - International Astronomical Union
//   IERS    - International Earth Rotation and Reference Systems Service
//
// Copyright 2026 AllMySat · Apache License 2.0
// ────────────────────────────────────────────────────────────────────────────

import Foundation

/// Fundamental physical and astronomical constants for astrodynamics.
///
/// Two reference systems are provided:
/// - **WGS-72**: Required by the SGP4 propagator (Spacetrack Report #3).
///   Using WGS-84 values with SGP4 introduces systematic errors.
/// - **WGS-84**: Modern geodetic standard, used for Keplerian calculations
///   outside the SGP4 context (e.g., TLE-derived orbital parameters).
///
public enum AstrodynamicsConstants {

    // MARK: - Earth - WGS-84

    /// Earth equatorial radius (WGS-84), in km.
    ///
    /// Source: NIMA TR8350.2, Table 3.1
    public static let earthRadiusKm: Double = 6378.137

    /// Earth mean radius (IUGG), in km.
    ///
    /// Arithmetic mean of equatorial and polar radii.
    /// Used for approximate altitude calculations from Keplerian elements.
    /// Source: Moritz, H., "Geodetic Reference System 1980", 2000
    public static let earthMeanRadiusKm: Double = 6371.0

    /// Earth gravitational parameter μ (WGS-84), in km³/s².
    ///
    /// μ = G × M_earth (product known to higher precision than G or M separately).
    /// Source: NIMA TR8350.2, Table 3.1
    public static let earthMu: Double = 398600.4418

    /// Earth rotation rate, in rad/s.
    ///
    /// Source: IERS Conventions (2010), Table 1.2
    public static let earthRotationRate: Double = 7.29211514670698e-5

    /// WGS-84 flattening (1/f).
    ///
    /// Source: NIMA TR8350.2, Table 3.1
    public static let earthFlattening: Double = 1.0 / 298.257223563

    // MARK: - Earth - WGS-72 (SGP4 only)

    /// Earth equatorial radius (WGS-72), in km.
    ///
    /// Used exclusively by the SGP4 propagator. Do not substitute WGS-84.
    /// Source: Hoots & Roehrich, Spacetrack Report #3, §3
    public static let earthRadiusWGS72: Double = 6378.135

    /// Earth gravitational parameter (WGS-72), in km³/s².
    ///
    /// Used exclusively by the SGP4 propagator. Do not substitute WGS-84.
    /// Source: Hoots & Roehrich, Spacetrack Report #3, §3
    public static let earthMuWGS72: Double = 398600.8

    /// WGS-72 flattening (1/f).
    ///
    /// Source: Hoots & Roehrich, Spacetrack Report #3
    public static let earthFlatteningWGS72: Double = 1.0 / 298.26

    // MARK: - J2000.0 Epoch

    /// Julian Date of the J2000.0 epoch (2000-01-01T12:00:00 UTC).
    ///
    /// Fundamental reference epoch for modern astrodynamics.
    /// Source: IAU 2009 System of Astronomical Constants
    public static let j2000JulianDate: Double = 2451545.0

    /// Unix timestamp of J2000.0 epoch, in seconds since 1970-01-01T00:00:00 UTC.
    public static let j2000UnixTimestamp: TimeInterval = 946728000.0

    // MARK: - Angle Conversions

    /// Degrees to radians.
    public static let deg2rad: Double = .pi / 180.0

    /// Radians to degrees.
    public static let rad2deg: Double = 180.0 / .pi

    // MARK: - Time

    /// Seconds in one solar day.
    public static let secondsPerDay: Double = 86400.0

    /// Minutes in one solar day.
    public static let minutesPerDay: Double = 1440.0

    // MARK: - Julian Date Conversion

    /// Converts a Swift `Date` to Julian Date.
    ///
    /// Uses the Meeus algorithm (Astronomical Algorithms, 2nd ed., Chapter 7).
    /// Valid for all dates in the Gregorian calendar.
    ///
    /// - Parameter date: The date to convert.
    /// - Returns: Julian Date as a fractional day count.
    public static func julianDate(from date: Date) -> Double {
        // Delegate to the SGP4Propagator's high-fidelity implementation
        // which uses Meeus calendar decomposition with Gregorian correction.
        SGP4Propagator.julianDate(from: date)
    }
}
