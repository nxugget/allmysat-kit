// ────────────────────────────────────────────────────────────────────────────
// TLEParser.swift
// AllMySatKit · TLE
//
// Parser for the NORAD Two-Line Element set format, extracting both raw
// orbital parameters and derived Keplerian quantities.
//
// A TLE consists of two 69-character lines, column-indexed:
//
//   Line 1:
//     Col  1      - Line number (1)
//     Col  3–7    - NORAD catalog number
//     Col  8      - Classification (U/C/S)
//     Col 10–17   - International designator (launch year + piece)
//     Col 19–32   - Epoch (2-digit year + fractional day of year)
//     Col 34–43   - First derivative of mean motion ÷ 2 (rev/day²)
//     Col 45–52   - Second derivative of mean motion ÷ 6 (implied decimal)
//     Col 54–61   - B* drag term (implied decimal + exponent)
//     Col 63      - Ephemeris type (0 = SGP4)
//     Col 65–68   - Element set number
//     Col 69      - Checksum (modulo 10)
//
//   Line 2:
//     Col  1      - Line number (2)
//     Col  3–7    - NORAD catalog number
//     Col  9–16   - Inclination (degrees)
//     Col 18–25   - RAAN (degrees)
//     Col 27–33   - Eccentricity (implied leading decimal point)
//     Col 35–42   - Argument of perigee (degrees)
//     Col 44–51   - Mean anomaly (degrees)
//     Col 53–63   - Mean motion (rev/day)
//     Col 64–68   - Revolution number at epoch
//     Col 69      - Checksum (modulo 10)
//
// Copyright 2026 AllMySat · Apache License 2.0
// ────────────────────────────────────────────────────────────────────────────

import Foundation

// MARK: - Orbital Elements

/// Keplerian orbital elements derived from a TLE.
///
/// Includes both the raw TLE fields and derived physical parameters
/// (semi-major axis, period, perigee/apogee altitudes).
public struct OrbitalElements: Sendable {

    // ── Raw TLE fields ─────────────────────────────────────

    /// NORAD catalog number (5-digit unique identifier).
    public let noradId: Int

    /// Classification: U (unclassified), C (classified), S (secret).
    public let classification: Character

    /// International designator (launch year, number, piece).
    public let internationalDesignator: String

    /// TLE epoch as a Julian Date.
    public let epochJulianDate: Double

    /// TLE epoch as a Swift `Date`.
    public let epochDate: Date

    /// Inclination, in degrees (0–180).
    public let inclination: Double

    /// Right Ascension of the Ascending Node, in degrees (0–360).
    public let raan: Double

    /// Orbital eccentricity (0 = circular, <1 = elliptical).
    public let eccentricity: Double

    /// Argument of perigee, in degrees (0–360).
    public let argumentOfPerigee: Double

    /// Mean anomaly at epoch, in degrees (0–360).
    public let meanAnomaly: Double

    /// Mean motion, in revolutions per day.
    public let meanMotion: Double

    /// B* drag coefficient (units: 1/Earth-radii).
    public let bstar: Double

    /// Revolution number at epoch.
    public let revolutionNumber: Int

    // ── Derived quantities ──────────────────────────────────

    /// Orbital period, in minutes.
    ///
    /// T = 1440 / n, where n = mean motion (rev/day)
    public let periodMinutes: Double

    /// Semi-major axis, in km.
    ///
    /// Derived from Kepler's Third Law:
    ///   a = (μ / n²)^(1/3)
    /// where n is in rad/s and μ = 398600.4418 km³/s² (WGS-84).
    public let semiMajorAxisKm: Double

    /// Perigee altitude above mean Earth radius, in km.
    ///
    ///   h_p = a(1 − e) − R_earth
    public let perigeeAltitudeKm: Double

    /// Apogee altitude above mean Earth radius, in km.
    ///
    ///   h_a = a(1 + e) − R_earth
    public let apogeeAltitudeKm: Double
}

// MARK: - TLE Parser

/// Parses NORAD Two-Line Element sets into structured orbital elements.
///
/// All methods are pure functions with no shared state.
public enum TLEParser {

    // MARK: - Public API

    /// Parses a TLE into orbital elements with derived quantities.
    ///
    /// - Parameters:
    ///   - line1: First line of the TLE (69 characters).
    ///   - line2: Second line of the TLE (69 characters).
    /// - Returns: Parsed orbital elements, or `nil` if the TLE is malformed.
    public static func parse(line1: String, line2: String) -> OrbitalElements? {
        guard line1.count >= 69, line2.count >= 69 else { return nil }

        let l1 = Array(line1)
        let l2 = Array(line2)

        // ── Line 1 fields ──────────────────────────────────────

        guard l1[0] == "1", l2[0] == "2" else { return nil }

        let noradId = Int(String(l1[2...6]).trimmingCharacters(in: .whitespaces)) ?? 0
        let classification = l1[7]
        let intlDesignator = String(l1[9...16]).trimmingCharacters(in: .whitespaces)

        // Epoch: 2-digit year + fractional day
        guard let epochYear = Int(String(l1[18...19]).trimmingCharacters(in: .whitespaces)),
              let epochDay = Double(String(l1[20...31]).trimmingCharacters(in: .whitespaces))
        else { return nil }

        let fullYear = epochYear < 57 ? 2000 + epochYear : 1900 + epochYear
        let epochJD = julianDateFromTLEEpoch(year: fullYear, day: epochDay)
        let epochDate = dateFromJulianDate(epochJD)

        // B* drag term: implied decimal + exponent
        let bstarStr = String(l1[53...60]).trimmingCharacters(in: .whitespaces)
        let bstar = parseImpliedDecimal(bstarStr)

        // ── Line 2 fields ──────────────────────────────────────

        guard let noradId2 = Int(String(l2[2...6]).trimmingCharacters(in: .whitespaces)),
              noradId2 == noradId else { return nil }

        guard let inclination = Double(String(l2[8...15]).trimmingCharacters(in: .whitespaces)),
              let raan = Double(String(l2[17...24]).trimmingCharacters(in: .whitespaces)),
              let argPerigee = Double(String(l2[34...41]).trimmingCharacters(in: .whitespaces)),
              let meanAnomaly = Double(String(l2[43...50]).trimmingCharacters(in: .whitespaces)),
              let meanMotion = Double(String(l2[52...62]).trimmingCharacters(in: .whitespaces))
        else { return nil }

        // Eccentricity: implied leading "0."
        let eccStr = "0." + String(l2[26...32]).trimmingCharacters(in: .whitespaces)
        guard let eccentricity = Double(eccStr) else { return nil }

        let revNumber = Int(String(l2[63...67]).trimmingCharacters(in: .whitespaces)) ?? 0

        // ── Derived quantities ─────────────────────────────────

        // Period: T = 1440 / n (minutes/day ÷ rev/day = minutes/rev)
        let periodMinutes = meanMotion > 0
            ? AstrodynamicsConstants.minutesPerDay / meanMotion
            : 0

        // Semi-major axis from Kepler's Third Law:
        //   n (rad/s) = meanMotion × 2π / 86400
        //   a = (μ / n²)^(1/3)
        let nRadPerSec = meanMotion * 2.0 * .pi / AstrodynamicsConstants.secondsPerDay
        let semiMajorAxis: Double
        if nRadPerSec > 0 {
            semiMajorAxis = pow(AstrodynamicsConstants.earthMu / (nRadPerSec * nRadPerSec),
                                1.0 / 3.0)
        } else {
            semiMajorAxis = 0
        }

        // Perigee and apogee altitudes
        let perigeeAlt = semiMajorAxis * (1.0 - eccentricity)
                        - AstrodynamicsConstants.earthMeanRadiusKm
        let apogeeAlt = semiMajorAxis * (1.0 + eccentricity)
                       - AstrodynamicsConstants.earthMeanRadiusKm

        return OrbitalElements(
            noradId: noradId,
            classification: classification,
            internationalDesignator: intlDesignator,
            epochJulianDate: epochJD,
            epochDate: epochDate,
            inclination: inclination,
            raan: raan,
            eccentricity: eccentricity,
            argumentOfPerigee: argPerigee,
            meanAnomaly: meanAnomaly,
            meanMotion: meanMotion,
            bstar: bstar,
            revolutionNumber: revNumber,
            periodMinutes: periodMinutes,
            semiMajorAxisKm: semiMajorAxis,
            perigeeAltitudeKm: perigeeAlt,
            apogeeAltitudeKm: apogeeAlt
        )
    }

    /// Validates a TLE line checksum.
    ///
    /// The checksum is the modulo-10 sum of all digits on the line,
    /// with '-' (minus signs) counted as 1.
    ///
    /// - Parameter line: A single TLE line (69 characters).
    /// - Returns: `true` if the checksum is valid.
    public static func validateChecksum(line: String) -> Bool {
        guard line.count >= 69 else { return false }

        let chars = Array(line)
        guard let expected = chars[68].wholeNumberValue else { return false }

        var sum = 0
        for i in 0..<68 {
            if let digit = chars[i].wholeNumberValue {
                sum += digit
            } else if chars[i] == "-" {
                sum += 1
            }
            // Letters and spaces contribute 0
        }

        return (sum % 10) == expected
    }

    /// Determines the age of a TLE relative to a reference date.
    ///
    /// - Parameters:
    ///   - elements: Parsed orbital elements containing the epoch.
    ///   - referenceDate: Date to measure age from (defaults to now).
    /// - Returns: Age in fractional days (positive = TLE is older than reference).
    public static func tleAge(elements: OrbitalElements, referenceDate: Date = Date()) -> Double {
        let refJD = AstrodynamicsConstants.julianDate(from: referenceDate)
        return refJD - elements.epochJulianDate
    }

    // MARK: - Internal Parsing

    /// Parses the implied-decimal exponential notation used for B* and nddot6.
    ///
    /// Format: " SMMMMM-N" where S = optional sign, M = mantissa digits,
    /// N = exponent. The decimal point is implied before the first digit.
    ///
    /// Example: " 11606-4" → 0.11606 × 10⁻⁴ = 0.0000011606
    private static func parseImpliedDecimal(_ field: String) -> Double {
        var s = field.trimmingCharacters(in: .whitespaces)
        guard !s.isEmpty else { return 0 }

        // Handle sign
        var sign: Double = 1.0
        if s.hasPrefix("-") {
            sign = -1.0
            s = String(s.dropFirst())
        } else if s.hasPrefix("+") {
            s = String(s.dropFirst())
        }

        // Split at the exponent sign
        let parts: [String]
        if let range = s.range(of: "-", range: s.index(after: s.startIndex)..<s.endIndex) {
            parts = [String(s[s.startIndex..<range.lowerBound]),
                     "-" + String(s[range.upperBound...])]
        } else if let range = s.range(of: "+", range: s.index(after: s.startIndex)..<s.endIndex) {
            parts = [String(s[s.startIndex..<range.lowerBound]),
                     String(s[range.upperBound...])]
        } else {
            return sign * (Double("0." + s) ?? 0)
        }

        let mantissa = Double("0." + parts[0]) ?? 0
        let exponent = Double(parts[1]) ?? 0
        return sign * mantissa * pow(10.0, exponent)
    }

    /// Converts a TLE epoch (year + fractional day) to a Julian Date.
    ///
    /// - Parameters:
    ///   - year: Full calendar year (e.g., 2024).
    ///   - day: Fractional day of year (1.0 = Jan 1 00:00 UTC).
    /// - Returns: Julian Date.
    private static func julianDateFromTLEEpoch(year: Int, day: Double) -> Double {
        // Build the epoch Date from TLE year + fractional day,
        // then delegate to the canonical Julian Date algorithm.
        var cal = Calendar(identifier: .gregorian)
        cal.timeZone = TimeZone(identifier: "UTC")!
        let jan1 = cal.date(from: DateComponents(year: year, month: 1, day: 1))!
        let epochDate = jan1.addingTimeInterval(
            (day - 1.0) * AstrodynamicsConstants.secondsPerDay
        )
        return AstrodynamicsConstants.julianDate(from: epochDate)
    }

    /// Converts a Julian Date back to a Swift `Date`.
    private static func dateFromJulianDate(_ jd: Double) -> Date {
        let unixTime = (jd - 2440587.5) * AstrodynamicsConstants.secondsPerDay
        return Date(timeIntervalSince1970: unixTime)
    }
}
