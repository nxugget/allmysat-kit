// ────────────────────────────────────────────────────────────────────────────
// PassPredictor.swift
// AllMySatKit · Propagation
//
// Satellite pass prediction over a ground observer.
//
// Computes Acquisition of Signal (AOS), Loss of Signal (LOS), and Time of
// Closest Approach (TCA) for satellite passes using a coarse-fine sweep
// with bisection-based horizon refinement.
//
// Algorithm:
//   1. Coarse scan (30 s steps) over the prediction window to detect
//      elevation sign changes (horizon crossings).
//   2. Fine scan (1 s steps) during active passes to sample the trajectory.
//   3. Binary search (15 iterations, ~0.001 s precision) at each horizon
//      crossing to pinpoint AOS and LOS times.
//   4. Filter passes by minimum peak elevation.
//
// The predictor is stateless and side-effect-free.
// Provide your own SGP4Propagator instance.
//
// Copyright 2026 AllMySat · Apache License 2.0
// ────────────────────────────────────────────────────────────────────────────

import Foundation
import CoreLocation

// MARK: - Result Types

/// A single point along a satellite's trajectory during a pass.
public struct TrajectoryPoint: Sendable {
    /// UTC time of this sample.
    public let time: Date
    /// Azimuth from observer, in degrees (0° = North, clockwise).
    public let azimuth: Double
    /// Elevation above observer's horizon, in degrees.
    public let elevation: Double
    /// Slant range from observer to satellite, in km.
    public let range: Double
    /// Range rate (positive = receding), in km/s.
    public let rangeRate: Double

    public init(time: Date, azimuth: Double, elevation: Double,
                range: Double, rangeRate: Double) {
        self.time = time
        self.azimuth = azimuth
        self.elevation = elevation
        self.range = range
        self.rangeRate = rangeRate
    }
}

/// A predicted satellite pass over a ground observer.
public struct PredictedPass: Sendable {
    /// Acquisition of Signal - moment the satellite rises above the horizon.
    public let aos: TrajectoryPoint
    /// Loss of Signal - moment the satellite sets below the horizon.
    public let los: TrajectoryPoint
    /// Time of Closest Approach - point of maximum elevation during the pass.
    public let tca: TrajectoryPoint
    /// Maximum elevation reached during the pass, in degrees.
    public let maxElevation: Double
    /// Sampled trajectory points (typically 1 s apart during the pass).
    public let trajectory: [TrajectoryPoint]

    /// Duration of the pass in seconds.
    public var duration: TimeInterval {
        los.time.timeIntervalSince(aos.time)
    }
}

// MARK: - Pass Predictor

/// Predicts satellite passes over a ground observer.
///
/// This is a pure computation module - no caching, no shared state, no singletons.
/// Provide an initialized `SGP4Propagator` and an observer location.
///
/// ## Usage
/// ```swift
/// let propagator = SGP4Propagator(tle1: line1, tle2: line2)!
/// let observer   = CLLocation(latitude: 48.8566, longitude: 2.3522)
/// let passes     = PassPredictor.predictPasses(
///     propagator: propagator,
///     observer: observer,
///     from: Date(),
///     to: Date().addingTimeInterval(86400),
///     minElevation: 10.0
/// )
/// ```
public enum PassPredictor {

    // MARK: - Configuration

    /// Coarse scan step (seconds). Chosen to be well below the minimum pass
    /// duration (~2 min for a 200 km LEO satellite) to avoid missing any pass.
    private static let coarseStep: TimeInterval = 30.0

    /// Fine scan step (seconds) during an active pass.
    private static let fineStep: TimeInterval = 1.0

    /// Number of bisection iterations for horizon refinement.
    /// 15 iterations: 30 s / 2¹⁵ ≈ 0.0009 s precision.
    private static let bisectionIterations = 15

    // MARK: - Public API

    /// Predicts all satellite passes within a time window.
    ///
    /// - Parameters:
    ///   - propagator: An initialized SGP4 propagator for the satellite.
    ///   - observer: Ground observer location (latitude, longitude, altitude).
    ///   - from: Start of the prediction window (UTC).
    ///   - to: End of the prediction window (UTC).
    ///   - minElevation: Minimum peak elevation to include a pass (degrees).
    ///     Passes whose maximum elevation is below this threshold are discarded.
    /// - Returns: Array of predicted passes, sorted chronologically.
    public static func predictPasses(
        propagator: SGP4Propagator,
        observer: CLLocation,
        from startDate: Date,
        to endDate: Date,
        minElevation: Double
    ) -> [PredictedPass] {

        var passes: [PredictedPass] = []
        var currentDate = startDate

        // State for the pass currently being tracked
        var inPass = false
        var passPoints: [TrajectoryPoint] = []
        var aosPoint: TrajectoryPoint?
        var tcaPoint: TrajectoryPoint?
        var maxEl: Double = 0

        while currentDate < endDate {
            guard let look = computeLookAngles(
                propagator: propagator, observer: observer, date: currentDate
            ) else {
                currentDate = currentDate.addingTimeInterval(coarseStep)
                continue
            }

            if look.elevation > 0 {
                // ── Inside a pass ──────────────────────────────────

                if !inPass {
                    inPass = true

                    // Refine the exact AOS moment via bisection
                    let refinedAOS = refineHorizonCrossing(
                        propagator: propagator, observer: observer,
                        before: currentDate.addingTimeInterval(-coarseStep),
                        after: currentDate,
                        rising: true
                    )
                    aosPoint = refinedAOS
                    passPoints = [refinedAOS]
                    maxEl = look.elevation
                    tcaPoint = TrajectoryPoint(
                        time: currentDate,
                        azimuth: look.azimuth,
                        elevation: look.elevation,
                        range: look.range,
                        rangeRate: look.rangeRate
                    )
                }

                let point = TrajectoryPoint(
                    time: currentDate,
                    azimuth: look.azimuth,
                    elevation: look.elevation,
                    range: look.range,
                    rangeRate: look.rangeRate
                )
                passPoints.append(point)

                if look.elevation > maxEl {
                    maxEl = look.elevation
                    tcaPoint = point
                }

                currentDate = currentDate.addingTimeInterval(fineStep)

            } else {
                // ── Outside a pass ─────────────────────────────────

                if inPass, let aos = aosPoint, let tca = tcaPoint {
                    // Refine the exact LOS moment via bisection
                    let refinedLOS = refineHorizonCrossing(
                        propagator: propagator, observer: observer,
                        before: currentDate.addingTimeInterval(-fineStep),
                        after: currentDate,
                        rising: false
                    )
                    passPoints.append(refinedLOS)

                    if maxEl >= minElevation {
                        passes.append(PredictedPass(
                            aos: aos,
                            los: refinedLOS,
                            tca: tca,
                            maxElevation: maxEl,
                            trajectory: passPoints
                        ))
                    }

                    inPass = false
                    aosPoint = nil
                    tcaPoint = nil
                    maxEl = 0
                    passPoints = []
                }

                currentDate = currentDate.addingTimeInterval(coarseStep)
            }
        }

        // Handle a pass that extends beyond the prediction window
        if inPass, let aos = aosPoint, let tca = tcaPoint, maxEl >= minElevation {
            let lastPoint = passPoints.last ?? tca
            passes.append(PredictedPass(
                aos: aos,
                los: lastPoint,
                tca: tca,
                maxElevation: maxEl,
                trajectory: passPoints
            ))
        }

        return passes
    }

    // MARK: - Horizon Refinement (Bisection)

    /// Finds the precise moment the satellite crosses the observer's horizon.
    ///
    /// Uses the bisection method on the elevation function:
    ///   f(t) = elevation(t)
    /// where f(before) and f(after) have opposite signs.
    ///
    /// - Parameters:
    ///   - propagator: SGP4 propagator instance.
    ///   - observer: Ground observer location.
    ///   - before: Time when satellite is on one side of the horizon.
    ///   - after: Time when satellite is on the other side.
    ///   - rising: `true` for AOS (elevation goes negative → positive),
    ///             `false` for LOS (positive → negative).
    /// - Returns: Trajectory point at the refined crossing time.
    private static func refineHorizonCrossing(
        propagator: SGP4Propagator,
        observer: CLLocation,
        before: Date,
        after: Date,
        rising: Bool
    ) -> TrajectoryPoint {
        var lo = before
        var hi = after
        var best = TrajectoryPoint(
            time: lo.addingTimeInterval(hi.timeIntervalSince(lo) / 2),
            azimuth: 0, elevation: 0, range: 0, rangeRate: 0
        )

        for _ in 0..<bisectionIterations {
            let mid = lo.addingTimeInterval(hi.timeIntervalSince(lo) / 2)
            guard let look = computeLookAngles(
                propagator: propagator, observer: observer, date: mid
            ) else { break }

            best = TrajectoryPoint(
                time: mid,
                azimuth: look.azimuth,
                elevation: look.elevation,
                range: look.range,
                rangeRate: look.rangeRate
            )

            // Bisect: narrow the interval containing the sign change
            if look.elevation > 0 {
                if rising { hi = mid } else { lo = mid }
            } else {
                if rising { lo = mid } else { hi = mid }
            }
        }

        return best
    }

    // MARK: - Internal

    /// Computes look angles from observer to satellite at a given time.
    private static func computeLookAngles(
        propagator: SGP4Propagator,
        observer: CLLocation,
        date: Date
    ) -> SGP4LookAngles? {
        guard let eci = propagator.propagate(date: date) else { return nil }
        return SGP4Propagator.lookAngles(
            position: eci.position,
            velocity: eci.velocity,
            observer: observer,
            date: date
        )
    }
}
