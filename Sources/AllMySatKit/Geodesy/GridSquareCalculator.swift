// ────────────────────────────────────────────────────────────────────────────
// GridSquareCalculator.swift
// AllMySatKit · Geodesy
//
// Maidenhead Locator System implementation for encoding geographic positions
// into grid square identifiers, and decoding them back.
//
// The Maidenhead system is a hierarchical geocoding scheme widely used in
// amateur radio to describe station locations with variable precision:
//
//   Precision   Characters   Grid Size        Example
//   ─────────   ──────────   ─────────        ───────
//   Field       2 (AA–RR)    20° × 10°        JN
//   Square      4 (AA00–RR99)  2° × 1°        JN18
//   Subsquare   6 (AA00aa–RR99xx)  5' × 2.5'  JN18eu
//   Extended    8 (AA00aa00–RR99xx99)  30" × 15"  JN18eu49
//
// The coordinate origin is at 90°S, 180°W (anti-meridian, South Pole).
// Longitude spans 360° across 18 fields (A–R), latitude 180° across 18 fields.
//
// Copyright 2026 AllMySat · Apache License 2.0
// ────────────────────────────────────────────────────────────────────────────

import Foundation
import CoreLocation

// MARK: - Precision

/// Level of precision for Maidenhead grid square encoding.
public enum GridPrecision: Int, Sendable, CaseIterable, Comparable {
    /// 2-character field (20° × 10°).
    case field = 2
    /// 4-character square (2° × 1°).
    case square = 4
    /// 6-character subsquare (5' × 2.5').
    case subsquare = 6
    /// 8-character extended subsquare (30" × 15").
    case extended = 8

    /// Infers precision from a locator string length.
    public static func from(locator: String) -> GridPrecision {
        switch locator.count {
        case 0..<4:  return .field
        case 4..<6:  return .square
        case 6..<8:  return .subsquare
        default:     return .extended
        }
    }

    /// Human-readable description of this precision level.
    public var displayName: String {
        switch self {
        case .field:     return "Field (2-char)"
        case .square:    return "Square (4-char)"
        case .subsquare: return "Subsquare (6-char)"
        case .extended:  return "Extended (8-char)"
        }
    }

    public static func < (lhs: GridPrecision, rhs: GridPrecision) -> Bool {
        lhs.rawValue < rhs.rawValue
    }
}

// MARK: - Grid Lines (for map overlays)

/// A line segment for rendering grid overlays on a map.
public struct GridLine: Sendable {
    /// Start coordinate of the line.
    public let start: CLLocationCoordinate2D
    /// End coordinate of the line.
    public let end: CLLocationCoordinate2D
    /// Grid precision level this line belongs to.
    public let precision: GridPrecision

    public init(start: CLLocationCoordinate2D, end: CLLocationCoordinate2D,
                precision: GridPrecision) {
        self.start = start
        self.end = end
        self.precision = precision
    }
}

// MARK: - Grid Square Calculator

/// Encodes and decodes Maidenhead grid square locators.
///
/// Performance note: an internal `NSCache` memoizes recent `coordinateToGridSquare`
/// results for repeated lookups with the same coordinates. The cache is thread-safe
/// and automatically evicts under memory pressure.
public enum GridSquareCalculator {

    /// Memoization cache for coordinate-to-grid-square lookups.
    ///
    /// NSCache is thread-safe and auto-evicts under memory pressure.
    /// Key = "lat,lon,precision", Value = locator string.
    private static let cache = NSCache<NSString, NSString>()

    // MARK: - Encode

    /// Converts a geographic coordinate to a Maidenhead grid square locator.
    ///
    /// The encoding algorithm:
    ///   1. Shift origin to (−90, −180) → all values positive
    ///   2. Field:     divide by 20° (lon) / 10° (lat) → uppercase A–R
    ///   3. Square:    divide remainder by 2° / 1°     → digit 0–9
    ///   4. Subsquare: divide remainder by 5' / 2.5'   → lowercase a–x
    ///   5. Extended:  divide remainder by 30" / 15"    → digit 0–9
    ///
    /// - Parameters:
    ///   - latitude: Latitude in decimal degrees (−90 to +90).
    ///   - longitude: Longitude in decimal degrees (−180 to +180).
    ///   - precision: Desired grid square precision.
    ///   - useCache: Whether to use the memoization cache (default: `true`).
    /// - Returns: Maidenhead locator string (e.g., `"JN18eu"`).
    public static func coordinateToGridSquare(
        latitude: Double,
        longitude: Double,
        precision: GridPrecision = .subsquare,
        useCache: Bool = true
    ) -> String {

        // Check cache
        let cacheKey = "\(latitude),\(longitude),\(precision.rawValue)" as NSString
        if useCache, let cached = cache.object(forKey: cacheKey) {
            return cached as String
        }

        // Clamp to valid ranges
        let lat = max(-90, min(90, latitude))
        let lon = max(-180, min(180, longitude))

        // Shift origin to (0, 0) at 90°S 180°W
        var adjLon = lon + 180.0
        var adjLat = lat + 90.0

        var result = ""

        // ── Field (2 chars): 20° lon × 10° lat ────────────────
        let fieldLon = Int(adjLon / 20.0)
        let fieldLat = Int(adjLat / 10.0)
        result.append(Character(UnicodeScalar(65 + min(fieldLon, 17))!))  // A–R
        result.append(Character(UnicodeScalar(65 + min(fieldLat, 17))!))  // A–R

        guard precision.rawValue > 2 else {
            cacheIfNeeded(result, forKey: cacheKey, useCache: useCache)
            return result
        }

        adjLon -= Double(fieldLon) * 20.0
        adjLat -= Double(fieldLat) * 10.0

        // ── Square (2 digits): 2° lon × 1° lat ────────────────
        let squareLon = Int(adjLon / 2.0)
        let squareLat = Int(adjLat / 1.0)
        result.append(Character(UnicodeScalar(48 + min(squareLon, 9))!))   // 0–9
        result.append(Character(UnicodeScalar(48 + min(squareLat, 9))!))   // 0–9

        guard precision.rawValue > 4 else {
            cacheIfNeeded(result, forKey: cacheKey, useCache: useCache)
            return result
        }

        adjLon -= Double(squareLon) * 2.0
        adjLat -= Double(squareLat) * 1.0

        // ── Subsquare (2 lowercase): 5' lon × 2.5' lat ────────
        // 2° / 24 = 5' = 0.0833°  |  1° / 24 = 2.5' = 0.0417°
        let subLon = Int(adjLon / (2.0 / 24.0))
        let subLat = Int(adjLat / (1.0 / 24.0))
        result.append(Character(UnicodeScalar(97 + min(subLon, 23))!))     // a–x
        result.append(Character(UnicodeScalar(97 + min(subLat, 23))!))     // a–x

        guard precision.rawValue > 6 else {
            cacheIfNeeded(result, forKey: cacheKey, useCache: useCache)
            return result
        }

        adjLon -= Double(subLon) * (2.0 / 24.0)
        adjLat -= Double(subLat) * (1.0 / 24.0)

        // ── Extended (2 digits): 30" lon × 15" lat ────────────
        // (2°/24)/10 = 30" = 0.00833°  |  (1°/24)/10 = 15" = 0.00417°
        let extLon = Int(adjLon / (2.0 / 240.0))
        let extLat = Int(adjLat / (1.0 / 240.0))
        result.append(Character(UnicodeScalar(48 + min(extLon, 9))!))      // 0–9
        result.append(Character(UnicodeScalar(48 + min(extLat, 9))!))      // 0–9

        cacheIfNeeded(result, forKey: cacheKey, useCache: useCache)
        return result
    }

    /// Batch conversion of coordinates to grid square locators.
    ///
    /// Bypasses the cache for efficiency when processing large datasets.
    ///
    /// - Parameters:
    ///   - coordinates: Array of (latitude, longitude) pairs.
    ///   - precision: Desired grid square precision.
    /// - Returns: Array of locator strings, in the same order as input.
    public static func batchCoordinateToGridSquare(
        coordinates: [(latitude: Double, longitude: Double)],
        precision: GridPrecision = .subsquare
    ) -> [String] {
        coordinates.map {
            coordinateToGridSquare(latitude: $0.latitude, longitude: $0.longitude,
                                   precision: precision, useCache: false)
        }
    }

    // MARK: - Decode

    /// Converts a Maidenhead locator back to the center coordinate of the grid cell.
    ///
    /// - Parameter locator: A valid Maidenhead locator string (2–8 characters).
    /// - Returns: The center of the grid cell, or `nil` if the locator is invalid.
    public static func gridSquareToCoordinate(locator: String) -> CLLocationCoordinate2D? {
        guard isValid(locator: locator) else { return nil }

        let chars = Array(locator.uppercased().unicodeScalars)

        // ── Field ──────────────────────────────────────────────
        var lon = (Double(chars[0].value) - 65.0) * 20.0 - 180.0
        var lat = (Double(chars[1].value) - 65.0) * 10.0 - 90.0
        var lonStep = 20.0
        var latStep = 10.0

        // ── Square ─────────────────────────────────────────────
        if chars.count >= 4 {
            lon += (Double(chars[2].value) - 48.0) * 2.0
            lat += (Double(chars[3].value) - 48.0) * 1.0
            lonStep = 2.0
            latStep = 1.0
        }

        // ── Subsquare ──────────────────────────────────────────
        if chars.count >= 6 {
            lon += (Double(chars[4].value) - 65.0) * (2.0 / 24.0)
            lat += (Double(chars[5].value) - 65.0) * (1.0 / 24.0)
            lonStep = 2.0 / 24.0
            latStep = 1.0 / 24.0
        }

        // ── Extended ───────────────────────────────────────────
        if chars.count >= 8 {
            lon += (Double(chars[6].value) - 48.0) * (2.0 / 240.0)
            lat += (Double(chars[7].value) - 48.0) * (1.0 / 240.0)
            lonStep = 2.0 / 240.0
            latStep = 1.0 / 240.0
        }

        // Return the center of the grid cell
        return CLLocationCoordinate2D(
            latitude: lat + latStep / 2.0,
            longitude: lon + lonStep / 2.0
        )
    }

    /// Returns the bounding box of a Maidenhead grid cell.
    ///
    /// - Parameter locator: A valid Maidenhead locator string.
    /// - Returns: Tuple of (minLat, maxLat, minLon, maxLon), or `nil` if invalid.
    public static func gridSquareBounds(locator: String)
        -> (minLat: Double, maxLat: Double, minLon: Double, maxLon: Double)? {

        guard let center = gridSquareToCoordinate(locator: locator) else { return nil }

        let precision = GridPrecision(rawValue: locator.count) ?? .subsquare
        let latStep: Double
        let lonStep: Double

        switch precision {
        case .field:    latStep = 10.0;          lonStep = 20.0
        case .square:   latStep = 1.0;           lonStep = 2.0
        case .subsquare: latStep = 1.0 / 24.0;  lonStep = 2.0 / 24.0
        case .extended: latStep = 1.0 / 240.0;  lonStep = 2.0 / 240.0
        }

        return (
            minLat: center.latitude - latStep / 2.0,
            maxLat: center.latitude + latStep / 2.0,
            minLon: center.longitude - lonStep / 2.0,
            maxLon: center.longitude + lonStep / 2.0
        )
    }

    // MARK: - Great-Circle Calculations

    /// Distance between two Maidenhead locators, in km.
    ///
    /// Uses `CLLocation.distance(from:)` which implements the WGS-84 ellipsoidal
    /// distance (Vincenty's formulae).
    ///
    /// - Returns: Distance in km, or `nil` if either locator is invalid.
    public static func distance(from locator1: String, to locator2: String) -> Double? {
        guard let coord1 = gridSquareToCoordinate(locator: locator1),
              let coord2 = gridSquareToCoordinate(locator: locator2) else { return nil }

        let loc1 = CLLocation(latitude: coord1.latitude, longitude: coord1.longitude)
        let loc2 = CLLocation(latitude: coord2.latitude, longitude: coord2.longitude)

        return loc1.distance(from: loc2) / 1000.0
    }

    /// Initial great-circle bearing from one locator to another, in degrees.
    ///
    /// Uses the forward azimuth formula:
    ///   θ = atan2( sin(Δλ)cos(φ₂), cos(φ₁)sin(φ₂) − sin(φ₁)cos(φ₂)cos(Δλ) )
    ///
    /// Source: Vincenty, T. - "Direct and Inverse Solutions of Geodesics on
    ///         the Ellipsoid with Application of Nested Equations", 1975
    ///         (simplified to spherical approximation)
    ///
    /// - Returns: Bearing in degrees (0° = North, clockwise), or `nil`.
    public static func bearing(from locator1: String, to locator2: String) -> Double? {
        guard let coord1 = gridSquareToCoordinate(locator: locator1),
              let coord2 = gridSquareToCoordinate(locator: locator2) else { return nil }

        let lat1 = coord1.latitude * AstrodynamicsConstants.deg2rad
        let lat2 = coord2.latitude * AstrodynamicsConstants.deg2rad
        let dLon = (coord2.longitude - coord1.longitude) * AstrodynamicsConstants.deg2rad

        let y = sin(dLon) * cos(lat2)
        let x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon)

        var bearing = atan2(y, x) * AstrodynamicsConstants.rad2deg
        if bearing < 0 { bearing += 360.0 }
        return bearing
    }

    // MARK: - Validation

    /// Validates a Maidenhead locator string.
    ///
    /// Accepted formats:
    ///   - 2 chars: field (AA–RR)
    ///   - 4 chars: square (AA00–RR99)
    ///   - 6 chars: subsquare (AA00aa–RR99xx)
    ///   - 8 chars: extended (AA00aa00–RR99xx99)
    ///
    /// - Parameter locator: String to validate.
    /// - Returns: `true` if the locator is syntactically valid.
    public static func isValid(locator: String) -> Bool {
        let pattern = "^[A-Ra-r]{2}([0-9]{2}([A-Xa-x]{2}([0-9]{2})?)?)?$"
        return locator.range(of: pattern, options: .regularExpression) != nil
    }

    // MARK: - Grid Lines for Map Overlay

    /// Generates grid lines for rendering Maidenhead boundaries on a map.
    ///
    /// - Parameters:
    ///   - minLat: Southern boundary of the visible region.
    ///   - maxLat: Northern boundary.
    ///   - minLon: Western boundary.
    ///   - maxLon: Eastern boundary.
    ///   - precision: Grid precision level to render.
    /// - Returns: Array of `GridLine` segments covering the visible region.
    public static func generateGridLines(
        minLat: Double, maxLat: Double,
        minLon: Double, maxLon: Double,
        precision: GridPrecision
    ) -> [GridLine] {

        let latStep: Double
        let lonStep: Double

        switch precision {
        case .field:     latStep = 10.0;          lonStep = 20.0
        case .square:    latStep = 1.0;           lonStep = 2.0
        case .subsquare: latStep = 1.0 / 24.0;   lonStep = 2.0 / 24.0
        case .extended:  latStep = 1.0 / 240.0;   lonStep = 2.0 / 240.0
        }

        var lines: [GridLine] = []

        // Snap to grid boundaries
        let startLat = floor(minLat / latStep) * latStep
        let startLon = floor(minLon / lonStep) * lonStep

        // Horizontal lines (constant latitude)
        var lat = startLat
        while lat <= maxLat {
            lines.append(GridLine(
                start: CLLocationCoordinate2D(latitude: lat, longitude: minLon),
                end: CLLocationCoordinate2D(latitude: lat, longitude: maxLon),
                precision: precision
            ))
            lat += latStep
        }

        // Vertical lines (constant longitude)
        var lon = startLon
        while lon <= maxLon {
            lines.append(GridLine(
                start: CLLocationCoordinate2D(latitude: minLat, longitude: lon),
                end: CLLocationCoordinate2D(latitude: maxLat, longitude: lon),
                precision: precision
            ))
            lon += lonStep
        }

        return lines
    }

    // MARK: - Internal

    private static func cacheIfNeeded(_ value: String, forKey key: NSString, useCache: Bool) {
        if useCache {
            cache.setObject(value as NSString, forKey: key)
        }
    }
}
