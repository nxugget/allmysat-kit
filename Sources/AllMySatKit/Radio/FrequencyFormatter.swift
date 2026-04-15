// ────────────────────────────────────────────────────────────────────────────
// FrequencyFormatter.swift
// AllMySatKit · Radio
//
// Handles radio frequency conversions, formatting, and ITU band
// classification for satellite communication frequencies.
//
// ITU Radio Frequency Band Designations (ITU-R V.431-8):
//   Band    Abbreviation    Range
//   ────    ────────────    ─────
//    4      VLF             3–30 kHz
//    5      LF              30–300 kHz
//    6      MF              300–3000 kHz
//    7      HF              3–30 MHz
//    8      VHF             30–300 MHz
//    9      UHF             300–3000 MHz
//   10      SHF             3–30 GHz
//   11      EHF             30–300 GHz
//
// Copyright 2026 AllMySat · Apache License 2.0
// ────────────────────────────────────────────────────────────────────────────

import Foundation

// MARK: - Display Unit

/// Frequency display unit for explicit formatting.
public enum FrequencyDisplayUnit: String, Sendable, CaseIterable {
    case hz  = "Hz"
    case khz = "kHz"
    case mhz = "MHz"
    case ghz = "GHz"

    /// Divisor to convert from Hz to this unit.
    var divisor: Double {
        switch self {
        case .hz:  return 1.0
        case .khz: return 1_000.0
        case .mhz: return 1_000_000.0
        case .ghz: return 1_000_000_000.0
        }
    }
}

// MARK: - Frequency Formatter

/// Formats and classifies radio frequencies.
///
/// All methods are stateless. Frequencies are represented as `Int64` Hz
/// internally for exact representation (no floating-point drift at GHz scale).
public enum FrequencyFormatter {

    // MARK: - Auto-Formatting

    /// Formats a frequency with automatic unit selection.
    ///
    /// Selects the most readable unit:
    ///   ≥ 1 GHz  → GHz (e.g., "10.450 GHz")
    ///   ≥ 1 MHz  → MHz (e.g., "145.800 MHz")
    ///   ≥ 1 kHz  → kHz (e.g., "67.000 kHz")
    ///   < 1 kHz  → Hz  (e.g., "500 Hz")
    ///
    /// - Parameter hz: Frequency in Hertz.
    /// - Returns: Human-readable frequency string with unit suffix.
    public static func format(hz: Int64) -> String {
        let absHz = abs(hz)
        if absHz >= 1_000_000_000 {
            return format(hz: hz, unit: .ghz)
        } else if absHz >= 1_000_000 {
            return format(hz: hz, unit: .mhz)
        } else if absHz >= 1_000 {
            return format(hz: hz, unit: .khz)
        } else {
            return format(hz: hz, unit: .hz)
        }
    }

    /// Formats a frequency in a specific unit.
    ///
    /// Uses 3 decimal places for MHz/GHz (sufficient for amateur radio
    /// channel spacing), 1 decimal place for kHz, and no decimals for Hz.
    ///
    /// - Parameters:
    ///   - hz: Frequency in Hertz.
    ///   - unit: Target display unit.
    /// - Returns: Formatted string with unit suffix.
    public static func format(hz: Int64, unit: FrequencyDisplayUnit) -> String {
        let value = Double(hz) / unit.divisor
        let format: String
        switch unit {
        case .ghz: format = "%.3f"
        case .mhz: format = "%.3f"
        case .khz: format = "%.1f"
        case .hz:  format = "%.0f"
        }
        return String(format: "\(format) %@", value, unit.rawValue)
    }

    // MARK: - SDR Export Format

    /// Formats a frequency for SDR (Software Defined Radio) software export.
    ///
    /// Outputs MHz with 6 decimal places and no unit suffix, matching the
    /// import format used by GQRX, SDR#, and SDR++.
    ///
    /// Example: 145800000 → "145.800000"
    ///
    /// - Parameter hz: Frequency in Hertz.
    /// - Returns: Frequency in MHz, 6 decimal places, no suffix.
    public static func formatForExport(hz: Int64) -> String {
        let mhz = Double(hz) / 1_000_000.0
        return String(format: "%.6f", mhz)
    }

    // MARK: - Parsing

    /// Parses a frequency string back to Hz.
    ///
    /// Accepts formats like:
    ///   "145.800 MHz", "10.450 GHz", "435000000", "67.000 kHz"
    ///
    /// The parser is case-insensitive for unit suffixes.
    ///
    /// - Parameter string: Frequency string to parse.
    /// - Returns: Frequency in Hz, or `nil` if the string is unparseable.
    public static func parse(_ string: String) -> Int64? {
        let trimmed = string.trimmingCharacters(in: .whitespaces)
        let lower = trimmed.lowercased()

        // Try to extract numeric value and unit
        let units: [(suffix: String, multiplier: Double)] = [
            ("ghz", 1_000_000_000),
            ("mhz", 1_000_000),
            ("khz", 1_000),
            ("hz", 1)
        ]

        for (suffix, multiplier) in units {
            if lower.hasSuffix(suffix) {
                let numStr = trimmed.dropLast(suffix.count)
                    .trimmingCharacters(in: .whitespaces)
                if let value = Double(numStr) {
                    return Int64(value * multiplier)
                }
            }
        }

        // Try bare numeric value (assume Hz)
        if let value = Double(trimmed) {
            return Int64(value)
        }

        return nil
    }

    // MARK: - ITU Band Classification

    /// Returns the ITU radio band name for a given frequency.
    ///
    /// Classification per ITU-R V.431-8:
    ///
    /// | Band | Abbreviation | Range           |
    /// |------|-------------|-----------------|
    /// |  4   | VLF         | 3–30 kHz        |
    /// |  5   | LF          | 30–300 kHz      |
    /// |  6   | MF          | 300–3,000 kHz   |
    /// |  7   | HF          | 3–30 MHz        |
    /// |  8   | VHF         | 30–300 MHz      |
    /// |  9   | UHF         | 300–3,000 MHz   |
    /// | 10   | SHF         | 3–30 GHz        |
    /// | 11   | EHF         | 30–300 GHz      |
    ///
    /// - Parameter hz: Frequency in Hertz.
    /// - Returns: ITU band abbreviation (e.g., "VHF"), or `nil` for out-of-range.
    public static func bandName(hz: Int64) -> String? {
        let absHz = abs(hz)
        switch absHz {
        case            3_000 ..<          30_000: return "VLF"
        case           30_000 ..<         300_000: return "LF"
        case          300_000 ..<       3_000_000: return "MF"
        case        3_000_000 ..<      30_000_000: return "HF"
        case       30_000_000 ..<     300_000_000: return "VHF"
        case      300_000_000 ..<   3_000_000_000: return "UHF"
        case    3_000_000_000 ..<  30_000_000_000: return "SHF"
        case   30_000_000_000 ..< 300_000_000_000: return "EHF"
        default: return nil
        }
    }

    // MARK: - Doppler Shift Formatting

    /// Formats a Doppler shift value with a sign prefix.
    ///
    /// Positive values indicate the satellite is receding (frequency shift down
    /// for the observer). Negative values indicate approach (shift up).
    ///
    /// - Parameter hz: Doppler shift in Hertz.
    /// - Returns: Formatted shift string (e.g., "+4.500 kHz", "−12.300 kHz").
    public static func formatShift(hz: Int64) -> String {
        let prefix = hz >= 0 ? "+" : ""
        let formatted = format(hz: hz)
        return "\(prefix)\(formatted)"
    }
}
