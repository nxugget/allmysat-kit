// ────────────────────────────────────────────────────────────────────────────
// CoordinateConverter.swift
// AllMySatKit · Geodesy
//
// Geodetic coordinate conversions: decimal degrees ↔ DMS (degrees, minutes,
// seconds), azimuth-elevation to unit vectors, compass directions, and
// distance unit conversions.
//
// DMS notation follows ISO 6709:
//   DD°MM'SS.ss"N    (latitude:  N/S suffix)
//   DDD°MM'SS.ss"E   (longitude: E/W suffix)
//
// The azimuth-to-direction conversion produces a unit vector in a
// right-handed coordinate system suitable for 3D rendering:
//   +X = East,  +Y = Up,  −Z = North
// This matches SceneKit/RealityKit conventions.
//
// Compass directions use the traditional 16-point wind rose (22.5° sectors)
// or the simplified 8-point rose (45° sectors).
//
// Copyright 2026 AllMySat · Apache License 2.0
// ────────────────────────────────────────────────────────────────────────────

import Foundation
import simd

// MARK: - Distance Unit

/// Unit of distance measurement.
public enum DistanceUnit: String, Sendable, CaseIterable {
    case kilometers
    case miles
    case nauticalMiles
}

// MARK: - Coordinate Converter

/// Pure geodetic coordinate conversions.
///
/// All methods are stateless, with no dependency on user preferences.
/// Formatting is the app's responsibility - this module provides only
/// the mathematical transformations.
public enum CoordinateConverter {

    // MARK: - Decimal ↔ DMS

    /// Converts a decimal degree value to degrees, minutes, seconds with
    /// a cardinal direction suffix.
    ///
    /// - Parameters:
    ///   - decimal: Coordinate value in decimal degrees.
    ///   - isLatitude: `true` for latitude (N/S), `false` for longitude (E/W).
    /// - Returns: Formatted DMS string, e.g. `"48°51'24.00"N"`.
    public static func toDMS(decimal: Double, isLatitude: Bool) -> String {
        let direction: String
        if isLatitude {
            direction = decimal >= 0 ? "N" : "S"
        } else {
            direction = decimal >= 0 ? "E" : "W"
        }

        let absolute = abs(decimal)
        let degrees = Int(absolute)
        let minutesFull = (absolute - Double(degrees)) * 60.0
        let minutes = Int(minutesFull)
        let seconds = (minutesFull - Double(minutes)) * 60.0

        return String(format: "%d°%02d'%05.2f\"%@", degrees, minutes, seconds, direction)
    }

    /// Converts latitude and longitude to a compact DMS string.
    ///
    /// - Returns: Combined string, e.g. `"48°51'24.00"N 2°21'07.00"E"`.
    public static func toDMSCompact(latitude: Double, longitude: Double) -> String {
        "\(toDMS(decimal: latitude, isLatitude: true)) \(toDMS(decimal: longitude, isLatitude: false))"
    }

    /// Converts degrees, minutes, seconds, and cardinal direction to decimal degrees.
    ///
    /// - Parameters:
    ///   - degrees: Whole degrees (non-negative).
    ///   - minutes: Arc-minutes (0–59).
    ///   - seconds: Arc-seconds (0–59.999...).
    ///   - direction: Cardinal direction: "N", "S", "E", or "W".
    /// - Returns: Decimal degrees. Southern and Western values are negative.
    public static func toDecimal(degrees: Int, minutes: Int, seconds: Double,
                                  direction: String) -> Double {
        let value = Double(degrees) + Double(minutes) / 60.0 + seconds / 3600.0
        let isNegative = direction.uppercased() == "S" || direction.uppercased() == "W"
        return isNegative ? -value : value
    }

    // MARK: - Azimuth/Elevation to Direction Vector

    /// Converts azimuth and elevation angles to a unit direction vector.
    ///
    /// Coordinate system (right-handed, matches SceneKit/RealityKit):
    ///   +X = East,  +Y = Up (zenith),  −Z = North
    ///
    /// The conversion from spherical to Cartesian:
    ///   x =  sin(az) × cos(el)
    ///   y =  sin(el)
    ///   z = −cos(az) × cos(el)
    ///
    /// - Parameters:
    ///   - azimuth: Azimuth in degrees (0° = North, clockwise).
    ///   - elevation: Elevation in degrees above the horizon.
    /// - Returns: Normalized 3D direction vector.
    public static func azElToDirection(azimuth: Double, elevation: Double) -> SIMD3<Float> {
        let azRad = azimuth * AstrodynamicsConstants.deg2rad
        let elRad = elevation * AstrodynamicsConstants.deg2rad

        let cosEl = cos(elRad)
        let x = sin(azRad) * cosEl
        let y = sin(elRad)
        let z = -cos(azRad) * cosEl

        return SIMD3<Float>(Float(x), Float(y), Float(z))
    }

    // MARK: - Distance Conversion

    /// Conversion factors from kilometers to other distance units.
    ///
    /// - km → miles:          × 0.621371 (exact: 1 / 1.609344)
    /// - km → nautical miles: × 0.539957 (exact: 1 / 1.852)
    ///
    /// Sources: BIPM SI Brochure, 9th ed., Table 8
    private static let conversionFactors: [DistanceUnit: Double] = [
        .kilometers: 1.0,
        .miles: 0.621371,
        .nauticalMiles: 0.539957
    ]

    /// Converts a distance from kilometers to the specified unit.
    ///
    /// - Parameters:
    ///   - km: Distance in kilometers.
    ///   - unit: Target distance unit.
    /// - Returns: Distance in the target unit.
    public static func convertDistance(_ km: Double, to unit: DistanceUnit) -> Double {
        km * (conversionFactors[unit] ?? 1.0)
    }

    // MARK: - Compass Direction (16-point)

    /// The 16-point compass rose, ordered from North clockwise.
    ///
    /// Each sector spans 22.5° centered on its cardinal/intercardinal direction.
    private static let compassPoints16 = [
        "N", "NNE", "NE", "ENE",
        "E", "ESE", "SE", "SSE",
        "S", "SSW", "SW", "WSW",
        "W", "WNW", "NW", "NNW"
    ]

    /// Returns the 16-point compass direction for a given azimuth.
    ///
    /// Uses 22.5° sectors with an 11.25° offset so that North spans
    /// [348.75°, 11.25°).
    ///
    /// - Parameter azimuth: Azimuth in degrees (0° = North, clockwise).
    /// - Returns: Compass direction abbreviation, e.g. `"NNE"`.
    public static func compassDirection(from azimuth: Double) -> String {
        // Normalize to [0, 360)
        var az = azimuth.truncatingRemainder(dividingBy: 360.0)
        if az < 0 { az += 360.0 }

        // 11.25° offset centers "N" on 0°, then divide by 22.5°
        let index = Int((az + 11.25).truncatingRemainder(dividingBy: 360.0) / 22.5)
        return compassPoints16[index % 16]
    }

    // MARK: - Compass Direction (8-point)

    /// The 8-point compass rose.
    private static let compassPoints8 = [
        "N", "NE", "E", "SE", "S", "SW", "W", "NW"
    ]

    /// Returns the 8-point compass direction for a given azimuth.
    ///
    /// Uses 45° sectors with a 22.5° offset.
    ///
    /// - Parameter azimuth: Azimuth in degrees (0° = North, clockwise).
    /// - Returns: Compass direction abbreviation, e.g. `"NE"`.
    public static func compassDirection8(from azimuth: Double) -> String {
        var az = azimuth.truncatingRemainder(dividingBy: 360.0)
        if az < 0 { az += 360.0 }

        let index = Int((az + 22.5).truncatingRemainder(dividingBy: 360.0) / 45.0)
        return compassPoints8[index % 8]
    }
}
