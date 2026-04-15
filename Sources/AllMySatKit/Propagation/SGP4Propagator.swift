// ────────────────────────────────────────────────────────────────────────────
// SGP4Propagator.swift
// AllMySatKit · Propagation
//
// Complete SGP4/SDP4 analytical orbital propagator for near-Earth and deep-space
// satellites based on Vallado et al. 2006 (AIAA 2006-6753).
//
// SGP4 handles near-Earth orbits (period < 225 minutes) with J₂, J₃, J₄
// perturbations and atmospheric drag. SDP4 handles deep-space satellites with
// Luni-Solar resonance coupling. Propagates NORAD Two-Line Element (TLE) sets
// to ECI state (position, velocity) and geodetic coordinates (lat, lon, altitude),
// with observer-relative look angles (azimuth, elevation, range, range-rate).
//
// Copyright 2026 AllMySat · Apache License 2.0
// ────────────────────────────────────────────────────────────────────────────

import Foundation
import CoreLocation
import simd

// MARK: - SGP4 Result Types

/// Earth-Centered Inertial state vector (TEME reference frame).
public struct SGP4ECIState: Sendable {
    public let position: SIMD3<Double>   // km
    public let velocity: SIMD3<Double>   // km/s

    public init(position: SIMD3<Double>, velocity: SIMD3<Double>) {
        self.position = position
        self.velocity = velocity
    }
}

/// Geodetic position on the Earth's surface.
public struct SGP4GeodeticPosition: Sendable {
    public let latitude: Double          // degrees  (-90 … +90)
    public let longitude: Double         // degrees  (-180 … +180)
    public let altitude: Double          // km above WGS-72 ellipsoid

    public var coordinate: CLLocationCoordinate2D {
        CLLocationCoordinate2D(latitude: latitude, longitude: longitude)
    }

    public init(latitude: Double, longitude: Double, altitude: Double) {
        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude
    }
}

/// Observer-relative look angles.
public struct SGP4LookAngles: Sendable {
    public let azimuth: Double           // degrees (0 = N, clockwise)
    public let elevation: Double         // degrees above horizon
    public let range: Double             // km
    public let rangeRate: Double         // km/s  (positive = receding)

    public init(azimuth: Double, elevation: Double, range: Double, rangeRate: Double) {
        self.azimuth = azimuth
        self.elevation = elevation
        self.range = range
        self.rangeRate = rangeRate
    }
}

// MARK: - SGP4Propagator

/// Full SGP4/SDP4 orbital propagator for near-Earth satellites.
///
/// Usage:
/// ```swift
/// guard let prop = SGP4Propagator(tle1: line1, tle2: line2) else { return }
/// guard let eci = prop.propagate(date: Date()) else { return }
/// let geo = SGP4Propagator.eciToGeodetic(position: eci.position, date: Date())
/// let look = SGP4Propagator.lookAngles(position: eci.position, velocity: eci.velocity,
///                                       observer: location, date: Date())
/// ```
/// 
/// WGS-72 constants used by SGP4 propagation, declared as nonisolated
/// to escape MainActor default isolation (these are pure constants).
nonisolated private enum C {
    static let twoPi       = 2.0 * Double.pi
    static let deg2rad     = Double.pi / 180.0
    static let rad2deg     = 180.0 / Double.pi
    static let minPerDay   = 1440.0
    static let secPerDay   = 86400.0
    static let x2o3        = 2.0 / 3.0

    static let earthRadius = 6378.135      // km  (WGS-72)
    static let mu          = 398600.8       // km³/s²  (WGS-72)
    static let xke         = 0.0743669161   // = 60/√(aE³/μ)  ER^1.5/min
    static let j2          = 0.001082616
    static let j3          = -0.00000253881
    static let j4          = -0.00000165597
    static let j3oj2       = j3 / j2
    static let vkmpersec   = earthRadius * xke / 60.0   // ER/min → km/s
    static let earthRotRate = 7.29211514670698e-5        // rad/s

    // Atmospheric model boundaries (km)
    static let qo: Double  = 120.0
    static let so: Double  = 78.0
}

public final class SGP4Propagator: @unchecked Sendable {

    // MARK: - Public Properties

    public let noradId: Int
    public let epochJd: Double
    public let inclinationDeg: Double
    public let periodMinutes: Double
    public let perigeeKm: Double
    public let isDeepSpace: Bool

    // MARK: - Parsed TLE Elements (radians / ER / rad·min⁻¹)

    private let bstar: Double
    nonisolated(unsafe) private var inclo: Double
    nonisolated(unsafe) private var nodeo: Double
    nonisolated(unsafe) private var ecco: Double
    nonisolated(unsafe) private var argpo: Double
    nonisolated(unsafe) private var mo: Double
    private let no_unkozai: Double      // Brouwer mean motion (rad/min)

    // MARK: - Initialization Constants (set once)

    // Trig
    private let cosio: Double
    private let cosio2: Double
    private let sinio: Double
    private let x1mth2: Double          // sin²i
    private let x7thm1: Double          // 7cos²i − 1
    private let con41: Double           // 3cos²i − 1

    // Drag coefficients
    private let eta: Double
    private let cc1: Double
    private let cc4: Double
    private let cc5: Double
    private let d2: Double
    private let d3: Double
    private let d4: Double

    // Secular rates
    private let mdot: Double            // rad/min
    private let argpdot: Double         // rad/min
    private let nodedot: Double         // rad/min
    private let nodecf: Double

    // Time coefficients
    private let t2cof: Double
    private let t3cof: Double
    private let t4cof: Double
    private let t5cof: Double

    // Periodic coefficients
    private let xlcof: Double
    private let aycof: Double
    private let omgcof: Double
    private let xmcof: Double
    private let delmo: Double
    private let sinmao: Double

    // Simplification flag (low-perigee drag model)
    private let isimp: Bool

    // MARK: - Deep-Space (SDP4) Stored Properties
    // Only populated for deep-space orbits (period ≥ 225 min).

    // Periodics coefficients (dscom → dpper)
    nonisolated(unsafe) private var ds_e3 = 0.0, ds_ee2 = 0.0
    nonisolated(unsafe) private var ds_peo = 0.0, ds_pgho = 0.0, ds_pho = 0.0
    nonisolated(unsafe) private var ds_pinco = 0.0, ds_plo = 0.0
    nonisolated(unsafe) private var ds_se2 = 0.0, ds_se3 = 0.0
    nonisolated(unsafe) private var ds_sgh2 = 0.0, ds_sgh3 = 0.0, ds_sgh4 = 0.0
    nonisolated(unsafe) private var ds_sh2 = 0.0, ds_sh3 = 0.0
    nonisolated(unsafe) private var ds_si2 = 0.0, ds_si3 = 0.0
    nonisolated(unsafe) private var ds_sl2 = 0.0, ds_sl3 = 0.0, ds_sl4 = 0.0
    nonisolated(unsafe) private var ds_xgh2 = 0.0, ds_xgh3 = 0.0, ds_xgh4 = 0.0
    nonisolated(unsafe) private var ds_xh2 = 0.0, ds_xh3 = 0.0
    nonisolated(unsafe) private var ds_xi2 = 0.0, ds_xi3 = 0.0
    nonisolated(unsafe) private var ds_xl2 = 0.0, ds_xl3 = 0.0, ds_xl4 = 0.0
    nonisolated(unsafe) private var ds_zmol = 0.0, ds_zmos = 0.0

    // Secular rates (dsinit → dspace)
    nonisolated(unsafe) private var ds_dedt = 0.0, ds_didt = 0.0, ds_dmdt = 0.0
    nonisolated(unsafe) private var ds_dnodt = 0.0, ds_domdt = 0.0

    // Resonance (dsinit → dspace)
    nonisolated(unsafe) private var ds_irez = 0
    nonisolated(unsafe) private var ds_d2201 = 0.0, ds_d2211 = 0.0
    nonisolated(unsafe) private var ds_d3210 = 0.0, ds_d3222 = 0.0
    nonisolated(unsafe) private var ds_d4410 = 0.0, ds_d4422 = 0.0
    nonisolated(unsafe) private var ds_d5220 = 0.0, ds_d5232 = 0.0
    nonisolated(unsafe) private var ds_d5421 = 0.0, ds_d5433 = 0.0
    nonisolated(unsafe) private var ds_del1 = 0.0, ds_del2 = 0.0, ds_del3 = 0.0
    nonisolated(unsafe) private var ds_xfact = 0.0, ds_xlamo = 0.0

    // Epoch sidereal time
    nonisolated(unsafe) private var ds_gsto = 0.0

    // Integration state (mutable during propagation)
    nonisolated(unsafe) private var ds_atime = 0.0, ds_xli = 0.0, ds_xni = 0.0

    // MARK: - Initialization

    /// Creates an SGP4 propagator from two TLE lines.
    /// Returns nil if the TLE cannot be parsed.
    /// Marked nonisolated since SGP4 is pure math with no main-thread dependencies.
    public nonisolated init?(tle1: String, tle2: String) {
        guard tle1.count >= 69, tle2.count >= 69 else { return nil }

        let l1 = Array(tle1)
        let l2 = Array(tle2)

        // ── Parse TLE ──────────────────────────────────────────────

        noradId = Int(String(l1[2...6]).trimmingCharacters(in: .whitespaces)) ?? 0

        // Epoch
        let epochYear = Int(String(l1[18...19]).trimmingCharacters(in: .whitespaces)) ?? 0
        let epochDay  = Double(String(l1[20...31]).trimmingCharacters(in: .whitespaces)) ?? 0
        let fullYear  = epochYear < 57 ? 2000 + epochYear : 1900 + epochYear
        epochJd = Self.julianDateFromEpoch(year: fullYear, day: epochDay)

        // B* drag
        let bstarStr = String(l1[53...60]).trimmingCharacters(in: .whitespaces)
        bstar = Self.parseExponentialField(bstarStr)

        // Line 2 orbital elements
        let incDeg  = Double(String(l2[8...15]).trimmingCharacters(in: .whitespaces))  ?? 0
        let raanDeg = Double(String(l2[17...24]).trimmingCharacters(in: .whitespaces)) ?? 0
        let eccStr  = "0." + String(l2[26...32]).trimmingCharacters(in: .whitespaces)
        let ecc     = Double(eccStr) ?? 0
        let argpDeg = Double(String(l2[34...41]).trimmingCharacters(in: .whitespaces)) ?? 0
        let maDeg   = Double(String(l2[43...50]).trimmingCharacters(in: .whitespaces)) ?? 0
        let mmRPD   = Double(String(l2[52...62]).trimmingCharacters(in: .whitespaces)) ?? 0

        inclo = incDeg  * C.deg2rad
        nodeo = raanDeg * C.deg2rad
        ecco  = ecc
        argpo = argpDeg * C.deg2rad
        mo    = maDeg   * C.deg2rad
        inclinationDeg = incDeg

        let no_kozai = mmRPD * C.twoPi / C.minPerDay   // rev/day → rad/min

        // ── Trig values ────────────────────────────────────────────

        cosio  = cos(inclo)
        cosio2 = cosio * cosio
        sinio  = sin(inclo)
        x1mth2 = 1.0 - cosio2
        x7thm1 = 7.0 * cosio2 - 1.0
        con41  = 3.0 * cosio2 - 1.0      // same as old x3thm1

        // ── Un-Kozai mean motion → Brouwer ─────────────────────────

        let eccsq  = ecco * ecco
        let omeosq = 1.0 - eccsq
        let rteosq = sqrt(omeosq)
        let cosio4 = cosio2 * cosio2

        let ak   = pow(C.xke / no_kozai, C.x2o3)
        let d1   = 0.75 * C.j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq)
        var del_  = d1 / (ak * ak)
        let adel  = ak * (1.0 - del_ * del_ - del_ * (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0))
        del_ = d1 / (adel * adel)
        no_unkozai = no_kozai / (1.0 + del_)
        let ao   = pow(C.xke / no_unkozai, C.x2o3)
        let posq = ao * ao

        // Perigee and period
        let rp = ao * (1.0 - ecco)
        perigeeKm     = (rp - 1.0) * C.earthRadius
        periodMinutes = C.twoPi / no_unkozai
        isDeepSpace   = periodMinutes >= 225.0

        // Simplification flag (low perigee or deep-space → simplified drag)
        isimp = rp < (220.0 / C.earthRadius + 1.0) || isDeepSpace

        // ── Atmospheric model parameters ───────────────────────────

        var sfour:  Double
        var qzms24: Double

        let qzms2ttemp = (C.qo - C.so) / C.earthRadius
        let qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp

        if perigeeKm < 156.0 {
            sfour = max(perigeeKm - 78.0, 20.0)
            let qztemp = (120.0 - sfour) / C.earthRadius
            qzms24 = qztemp * qztemp * qztemp * qztemp
            sfour = sfour / C.earthRadius + 1.0
        } else {
            sfour  = C.so / C.earthRadius + 1.0
            qzms24 = qzms2t
        }

        let pinvsq = 1.0 / posq
        let tsi    = 1.0 / (ao - sfour)
        eta = ao * ecco * tsi
        let etasq = eta * eta
        let eeta  = ecco * eta
        let psisq = abs(1.0 - etasq)
        let coef  = qzms24 * pow(tsi, 4.0)
        let coef1 = coef / pow(psisq, 3.5)

        // ── C1 – C5 ───────────────────────────────────────────────

        let cc2 = coef1 * no_unkozai *
            (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
             + 0.375 * C.j2 * tsi / psisq * con41
               * (8.0 + 3.0 * etasq * (8.0 + etasq)))

        cc1 = bstar * cc2

        var cc3_: Double = 0
        if ecco > 1.0e-4 {
            cc3_ = -2.0 * coef * tsi * C.j3oj2 * no_unkozai * sinio
        }

        cc4 = 2.0 * no_unkozai * coef1 * ao * omeosq *
            (eta * (2.0 + 0.5 * etasq) + ecco * (0.5 + 2.0 * etasq)
             - C.j2 * tsi / (ao * psisq)
               * (-3.0 * con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta))
                  + 0.75 * x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * argpo)))

        cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq)

        // ── Secular rates (J2, J4) ────────────────────────────────

        let temp1 = 1.5 * C.j2 * pinvsq * no_unkozai
        let temp2 = 0.5 * temp1 * C.j2 * pinvsq
        let temp3 = -0.46875 * C.j4 * pinvsq * pinvsq * no_unkozai

        mdot = no_unkozai
            + 0.5 * temp1 * rteosq * con41
            + 0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4)

        argpdot = -0.5 * temp1 * (1.0 - 5.0 * cosio2)     // NB: con42 = 1 - 5cos²i
            + 0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4)
            + temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4)

        let xhdot1 = -temp1 * cosio
        nodedot = xhdot1
            + (0.5 * temp2 * (4.0 - 19.0 * cosio2)
               + 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio

        // ── Drag terms ─────────────────────────────────────────────

        nodecf = 3.5 * omeosq * xhdot1 * cc1
        t2cof  = 1.5 * cc1

        omgcof = bstar * cc3_ * cos(argpo)

        if ecco > 1.0e-4 {
            xmcof = -C.x2o3 * coef * bstar / eeta
        } else {
            xmcof = 0
        }

        let delmotemp = 1.0 + eta * cos(mo)
        delmo  = delmotemp * delmotemp * delmotemp
        sinmao = sin(mo)

        // Long-period periodics
        let temp4 = 1.5e-12
        if abs(cosio + 1.0) > temp4 {
            xlcof = -0.25 * C.j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio)
        } else {
            xlcof = -0.25 * C.j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4
        }
        aycof = -0.5 * C.j3oj2 * sinio

        // x7thm1 already set above

        // ── Deep-space (SDP4) initialization ────────────────────

        if isDeepSpace {
            let xpidot = argpdot + nodedot
            ds_gsto = Self.gmst(jd: epochJd)
            let dsEpoch = epochJd - 2433281.5   // days from Jan 0 1950

            // ── dscom: solar/lunar coupling coefficients ────────

            let zes: Double  = 0.01675
            let zel: Double  = 0.05490
            let c1ss: Double = 2.9864797e-6
            let c1l: Double  = 4.7968065e-7
            let zsinis: Double = 0.39785416
            let zcosis: Double = 0.91744867
            let zcosgs: Double = 0.1945905
            let zsings: Double = -0.98088458

            var dc_nm   = no_unkozai
            var dc_em   = ecco
            let dc_snodm  = sin(nodeo)
            let dc_cnodm  = cos(nodeo)
            let dc_sinomm = sin(argpo)
            let dc_cosomm = cos(argpo)
            let dc_sinim  = sin(inclo)
            let dc_cosim  = cos(inclo)
            let dc_emsq   = dc_em * dc_em
            let dc_betasq = 1.0 - dc_emsq
            let dc_rtemsq = sqrt(dc_betasq)

            ds_peo = 0; ds_pinco = 0; ds_plo = 0; ds_pgho = 0; ds_pho = 0
            let dc_day = dsEpoch + 18261.5
            let dc_xnodce = (4.5236020 - 9.2422029e-4 * dc_day)
                .truncatingRemainder(dividingBy: C.twoPi)
            let dc_stem = sin(dc_xnodce)
            let dc_ctem = cos(dc_xnodce)
            let dc_zcosil = 0.91375164 - 0.03568096 * dc_ctem
            let dc_zsinil = sqrt(1.0 - dc_zcosil * dc_zcosil)
            let dc_zsinhl = 0.089683511 * dc_stem / dc_zsinil
            let dc_zcoshl = sqrt(1.0 - dc_zsinhl * dc_zsinhl)
            let dc_gam    = 5.8351514 + 0.0019443680 * dc_day
            var dc_zx     = 0.39785416 * dc_stem / dc_zsinil
            let dc_zy     = dc_zcoshl * dc_ctem + 0.91744867 * dc_zsinhl * dc_stem
            dc_zx = atan2(dc_zx, dc_zy)
            dc_zx = dc_gam + dc_zx - dc_xnodce
            let dc_zcosgl = cos(dc_zx)
            let dc_zsingl = sin(dc_zx)

            var dc_zcosg = zcosgs; var dc_zsing = zsings
            var dc_zcosi = zcosis; var dc_zsini = zsinis
            var dc_zcosh = dc_cnodm; var dc_zsinh = dc_snodm
            var dc_cc    = c1ss
            let dc_xnoi  = 1.0 / dc_nm

            var dc_ss1 = 0.0, dc_ss2 = 0.0, dc_ss3 = 0.0, dc_ss4 = 0.0
            var dc_ss5 = 0.0, dc_ss6 = 0.0, dc_ss7 = 0.0
            var dc_sz1 = 0.0, dc_sz2 = 0.0, dc_sz3 = 0.0
            var dc_sz11 = 0.0, dc_sz12 = 0.0, dc_sz13 = 0.0
            var dc_sz21 = 0.0, dc_sz22 = 0.0, dc_sz23 = 0.0
            var dc_sz31 = 0.0, dc_sz32 = 0.0, dc_sz33 = 0.0
            var dc_s1 = 0.0, dc_s2 = 0.0, dc_s3 = 0.0, dc_s4 = 0.0
            var dc_s5 = 0.0, dc_s6 = 0.0, dc_s7 = 0.0
            var dc_z1 = 0.0, dc_z2 = 0.0, dc_z3 = 0.0
            var dc_z11 = 0.0, dc_z12 = 0.0, dc_z13 = 0.0
            var dc_z21 = 0.0, dc_z22 = 0.0, dc_z23 = 0.0
            var dc_z31 = 0.0, dc_z32 = 0.0, dc_z33 = 0.0

            for lsflg in 1...2 {
                let a1  =  dc_zcosg * dc_zcosh + dc_zsing * dc_zcosi * dc_zsinh
                let a3  = -dc_zsing * dc_zcosh + dc_zcosg * dc_zcosi * dc_zsinh
                let a7  = -dc_zcosg * dc_zsinh + dc_zsing * dc_zcosi * dc_zcosh
                let a8  =  dc_zsing * dc_zsini
                let a9  =  dc_zsing * dc_zsinh + dc_zcosg * dc_zcosi * dc_zcosh
                let a10 =  dc_zcosg * dc_zsini
                let a2  =  dc_cosim * a7 + dc_sinim * a8
                let a4  =  dc_cosim * a9 + dc_sinim * a10
                let a5  = -dc_sinim * a7 + dc_cosim * a8
                let a6  = -dc_sinim * a9 + dc_cosim * a10

                let x1 =  a1 * dc_cosomm + a2 * dc_sinomm
                let x2 =  a3 * dc_cosomm + a4 * dc_sinomm
                let x3 = -a1 * dc_sinomm + a2 * dc_cosomm
                let x4 = -a3 * dc_sinomm + a4 * dc_cosomm
                let x5 =  a5 * dc_sinomm
                let x6 =  a6 * dc_sinomm
                let x7 =  a5 * dc_cosomm
                let x8 =  a6 * dc_cosomm

                dc_z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3
                dc_z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4
                dc_z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4
                dc_z1  = 3.0 * (a1 * a1 + a2 * a2) + dc_z31 * dc_emsq
                dc_z2  = 6.0 * (a1 * a3 + a2 * a4) + dc_z32 * dc_emsq
                dc_z3  = 3.0 * (a3 * a3 + a4 * a4) + dc_z33 * dc_emsq
                dc_z11 = -6.0 * a1 * a5 + dc_emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5)
                dc_z12 = -6.0 * (a1 * a6 + a3 * a5) + dc_emsq *
                    (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5))
                dc_z13 = -6.0 * a3 * a6 + dc_emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6)
                dc_z21 = 6.0 * a2 * a5 + dc_emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7)
                dc_z22 = 6.0 * (a4 * a5 + a2 * a6) + dc_emsq *
                    (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8))
                dc_z23 = 6.0 * a4 * a6 + dc_emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8)
                dc_z1  = dc_z1 + dc_z1 + dc_betasq * dc_z31
                dc_z2  = dc_z2 + dc_z2 + dc_betasq * dc_z32
                dc_z3  = dc_z3 + dc_z3 + dc_betasq * dc_z33
                dc_s3  = dc_cc * dc_xnoi
                dc_s2  = -0.5 * dc_s3 / dc_rtemsq
                dc_s4  = dc_s3 * dc_rtemsq
                dc_s1  = -15.0 * dc_em * dc_s4
                dc_s5  = x1 * x3 + x2 * x4
                dc_s6  = x2 * x3 + x1 * x4
                dc_s7  = x2 * x4 - x1 * x3

                if lsflg == 1 {
                    dc_ss1 = dc_s1; dc_ss2 = dc_s2; dc_ss3 = dc_s3; dc_ss4 = dc_s4
                    dc_ss5 = dc_s5; dc_ss6 = dc_s6; dc_ss7 = dc_s7
                    dc_sz1 = dc_z1; dc_sz2 = dc_z2; dc_sz3 = dc_z3
                    dc_sz11 = dc_z11; dc_sz12 = dc_z12; dc_sz13 = dc_z13
                    dc_sz21 = dc_z21; dc_sz22 = dc_z22; dc_sz23 = dc_z23
                    dc_sz31 = dc_z31; dc_sz32 = dc_z32; dc_sz33 = dc_z33
                    dc_zcosg = dc_zcosgl; dc_zsing = dc_zsingl
                    dc_zcosi = dc_zcosil; dc_zsini = dc_zsinil
                    dc_zcosh = dc_zcoshl * dc_cnodm + dc_zsinhl * dc_snodm
                    dc_zsinh = dc_snodm * dc_zcoshl - dc_cnodm * dc_zsinhl
                    dc_cc = c1l
                }
            }

            ds_zmol = (4.7199672 + 0.22997150 * dc_day - dc_gam)
                .truncatingRemainder(dividingBy: C.twoPi)
            ds_zmos = (6.2565837 + 0.017201977 * dc_day)
                .truncatingRemainder(dividingBy: C.twoPi)

            // Solar terms
            ds_se2  =  2.0 * dc_ss1 * dc_ss6
            ds_se3  =  2.0 * dc_ss1 * dc_ss7
            ds_si2  =  2.0 * dc_ss2 * dc_sz12
            ds_si3  =  2.0 * dc_ss2 * (dc_sz13 - dc_sz11)
            ds_sl2  = -2.0 * dc_ss3 * dc_sz2
            ds_sl3  = -2.0 * dc_ss3 * (dc_sz3 - dc_sz1)
            ds_sl4  = -2.0 * dc_ss3 * (-21.0 - 9.0 * dc_emsq) * zes
            ds_sgh2 =  2.0 * dc_ss4 * dc_sz32
            ds_sgh3 =  2.0 * dc_ss4 * (dc_sz33 - dc_sz31)
            ds_sgh4 = -18.0 * dc_ss4 * zes
            ds_sh2  = -2.0 * dc_ss2 * dc_sz22
            ds_sh3  = -2.0 * dc_ss2 * (dc_sz23 - dc_sz21)

            // Lunar terms
            ds_ee2  =  2.0 * dc_s1 * dc_s6
            ds_e3   =  2.0 * dc_s1 * dc_s7
            ds_xi2  =  2.0 * dc_s2 * dc_z12
            ds_xi3  =  2.0 * dc_s2 * (dc_z13 - dc_z11)
            ds_xl2  = -2.0 * dc_s3 * dc_z2
            ds_xl3  = -2.0 * dc_s3 * (dc_z3 - dc_z1)
            ds_xl4  = -2.0 * dc_s3 * (-21.0 - 9.0 * dc_emsq) * zel
            ds_xgh2 =  2.0 * dc_s4 * dc_z32
            ds_xgh3 =  2.0 * dc_s4 * (dc_z33 - dc_z31)
            ds_xgh4 = -18.0 * dc_s4 * zel
            ds_xh2  = -2.0 * dc_s2 * dc_z22
            ds_xh3  = -2.0 * dc_s2 * (dc_z23 - dc_z21)

            // ── dpper (init): adjust mean elements at epoch ─────

            var dp_zm = ds_zmos
            var dp_zf = dp_zm + 2.0 * zes * sin(dp_zm)
            var dp_sinzf = sin(dp_zf)
            var dp_f2 = 0.5 * dp_sinzf * dp_sinzf - 0.25
            var dp_f3 = -0.5 * dp_sinzf * cos(dp_zf)
            let dp_ses = ds_se2 * dp_f2 + ds_se3 * dp_f3
            let dp_sis = ds_si2 * dp_f2 + ds_si3 * dp_f3
            let dp_sls = ds_sl2 * dp_f2 + ds_sl3 * dp_f3 + ds_sl4 * dp_sinzf
            let dp_sghs = ds_sgh2 * dp_f2 + ds_sgh3 * dp_f3 + ds_sgh4 * dp_sinzf
            let dp_shs = ds_sh2 * dp_f2 + ds_sh3 * dp_f3

            dp_zm = ds_zmol
            dp_zf = dp_zm + 2.0 * zel * sin(dp_zm)
            dp_sinzf = sin(dp_zf)
            dp_f2 = 0.5 * dp_sinzf * dp_sinzf - 0.25
            dp_f3 = -0.5 * dp_sinzf * cos(dp_zf)
            let dp_sel  = ds_ee2 * dp_f2 + ds_e3 * dp_f3
            let dp_sil  = ds_xi2 * dp_f2 + ds_xi3 * dp_f3
            let dp_sll  = ds_xl2 * dp_f2 + ds_xl3 * dp_f3 + ds_xl4 * dp_sinzf
            let dp_sghl = ds_xgh2 * dp_f2 + ds_xgh3 * dp_f3 + ds_xgh4 * dp_sinzf
            let dp_shll = ds_xh2 * dp_f2 + ds_xh3 * dp_f3

            let dp_pe   = dp_ses + dp_sel
            let dp_pinc = dp_sis + dp_sil
            let dp_pl   = dp_sls + dp_sll
            var dp_pgh  = dp_sghs + dp_sghl
            var dp_ph   = dp_shs + dp_shll

            ds_peo   = dp_pe
            ds_pinco = dp_pinc
            ds_plo   = dp_pl
            ds_pgho  = dp_pgh
            ds_pho   = dp_ph

            inclo  = inclo + dp_pinc
            ecco   = ecco + dp_pe
            let dp_sinip = sin(inclo)
            let dp_cosip = cos(inclo)

            if inclo >= 0.2 {
                dp_ph  = dp_ph / dp_sinip
                dp_pgh = dp_pgh - dp_cosip * dp_ph
                argpo  = argpo + dp_pgh
                nodeo  = nodeo + dp_ph
                mo     = mo + dp_pl
            } else {
                let dp_sinop = sin(nodeo)
                let dp_cosop = cos(nodeo)
                var dp_alfdp = dp_sinip * dp_sinop
                var dp_betdp = dp_sinip * dp_cosop
                let dp_dalf  =  dp_ph * dp_cosop + dp_pinc * dp_cosip * dp_sinop
                let dp_dbet  = -dp_ph * dp_sinop + dp_pinc * dp_cosip * dp_cosop
                dp_alfdp = dp_alfdp + dp_dalf
                dp_betdp = dp_betdp + dp_dbet
                nodeo = nodeo.truncatingRemainder(dividingBy: C.twoPi)
                if nodeo < 0 { nodeo += C.twoPi }
                let dp_xls = mo + argpo + dp_pl + dp_pgh
                    + (dp_cosip - dp_pinc * dp_sinip) * nodeo
                let dp_xnoh = nodeo
                nodeo = atan2(dp_alfdp, dp_betdp)
                if abs(dp_xnoh - nodeo) > Double.pi {
                    nodeo = nodeo < dp_xnoh ? nodeo + C.twoPi : nodeo - C.twoPi
                }
                mo = mo + dp_pl
                argpo = dp_xls - mo - dp_cosip * nodeo
            }

            // ── dsinit: resonance and secular rates ─────────────

            let znl_di: Double = 1.5835218e-4
            let zns_di: Double = 1.19459e-5
            let rptim: Double  = 4.37526908801129966e-3

            // Solar secular effects
            let di_ses = dc_ss1 * zns_di * dc_ss5
            let di_sis = dc_ss2 * zns_di * (dc_sz11 + dc_sz13)
            let di_sls = -zns_di * dc_ss3 * (dc_sz1 + dc_sz3 - 14.0 - 6.0 * dc_emsq)
            let di_sghs = dc_ss4 * zns_di * (dc_sz31 + dc_sz33 - 6.0)
            var di_shs = -zns_di * dc_ss2 * (dc_sz21 + dc_sz23)
            if inclo < 5.2359877e-2 || inclo > Double.pi - 5.2359877e-2 {
                di_shs = 0
            }
            if dc_sinim != 0 { di_shs = di_shs / dc_sinim }
            let di_sgs = di_sghs - dc_cosim * di_shs

            // Lunar secular effects
            ds_dedt = di_ses + dc_s1 * znl_di * dc_s5
            ds_didt = di_sis + dc_s2 * znl_di * (dc_z11 + dc_z13)
            ds_dmdt = di_sls - znl_di * dc_s3 * (dc_z1 + dc_z3 - 14.0 - 6.0 * dc_emsq)
            let di_sghl = dc_s4 * znl_di * (dc_z31 + dc_z33 - 6.0)
            var di_shll = -znl_di * dc_s2 * (dc_z21 + dc_z23)
            if inclo < 5.2359877e-2 || inclo > Double.pi - 5.2359877e-2 {
                di_shll = 0
            }
            ds_domdt = di_sgs + di_sghl
            ds_dnodt = di_shs
            if dc_sinim != 0 {
                ds_domdt = ds_domdt - dc_cosim / dc_sinim * di_shll
                ds_dnodt = ds_dnodt + di_shll / dc_sinim
            }

            // Resonance flag
            ds_irez = 0
            if dc_nm >= 0.0034906585 && dc_nm < 0.0052359877 { ds_irez = 1 }
            if dc_nm >= 8.26e-3 && dc_nm <= 9.24e-3 && dc_em >= 0.5 { ds_irez = 2 }

            // ── Geopotential resonance for 12-hour orbits ───────

            if ds_irez == 2 {
                let di_cosisq = dc_cosim * dc_cosim
                let di_emo    = dc_em
                dc_em = ecco
                let di_emsq_  = ecco * ecco
                let di_eoc    = dc_em * di_emsq_
                let di_g201   = -0.306 - (dc_em - 0.64) * 0.440

                var di_g211: Double, di_g310: Double, di_g322: Double
                var di_g410: Double, di_g422: Double, di_g520: Double

                if dc_em <= 0.65 {
                    di_g211 =    3.616  -  13.2470 * dc_em +  16.2900 * di_emsq_
                    di_g310 =  -19.302  + 117.3900 * dc_em - 228.4190 * di_emsq_ + 156.5910 * di_eoc
                    di_g322 =  -18.9068 + 109.7927 * dc_em - 214.6334 * di_emsq_ + 146.5816 * di_eoc
                    di_g410 =  -41.122  + 242.6940 * dc_em - 471.0940 * di_emsq_ + 313.9530 * di_eoc
                    di_g422 = -146.407  + 841.8800 * dc_em - 1629.014 * di_emsq_ + 1083.4350 * di_eoc
                    di_g520 = -532.114  + 3017.977 * dc_em - 5740.032 * di_emsq_ + 3708.2760 * di_eoc
                } else {
                    di_g211 =   -72.099 +   331.819 * dc_em -   508.738 * di_emsq_ +   266.724 * di_eoc
                    di_g310 =  -346.844 +  1582.851 * dc_em -  2415.925 * di_emsq_ +  1246.113 * di_eoc
                    di_g322 =  -342.585 +  1554.908 * dc_em -  2366.899 * di_emsq_ +  1215.972 * di_eoc
                    di_g410 = -1052.797 +  4758.686 * dc_em -  7193.992 * di_emsq_ +  3651.957 * di_eoc
                    di_g422 = -3581.690 + 16178.110 * dc_em - 24462.770 * di_emsq_ + 12422.520 * di_eoc
                    di_g520 = dc_em > 0.715
                        ? -5149.66 + 29936.92 * dc_em - 54087.36 * di_emsq_ + 31324.56 * di_eoc
                        :  1464.74 -  4664.75 * dc_em +  3763.64 * di_emsq_
                }

                var di_g533: Double, di_g521: Double, di_g532: Double
                if dc_em < 0.7 {
                    di_g533 = -919.22770 + 4988.6100 * dc_em - 9064.7700 * di_emsq_ + 5542.21  * di_eoc
                    di_g521 = -822.71072 + 4568.6173 * dc_em - 8491.4146 * di_emsq_ + 5337.524 * di_eoc
                    di_g532 = -853.66600 + 4690.2500 * dc_em - 8624.7700 * di_emsq_ + 5341.4   * di_eoc
                } else {
                    di_g533 = -37995.780 + 161616.52  * dc_em - 229838.20  * di_emsq_ + 109377.94 * di_eoc
                    di_g521 = -51752.104 + 218913.95  * dc_em - 309468.16  * di_emsq_ + 146349.42 * di_eoc
                    di_g532 = -40023.880 + 170470.89  * dc_em - 242699.48  * di_emsq_ + 115605.82 * di_eoc
                }

                let di_sini2 = dc_sinim * dc_sinim
                let di_f220 = 0.75 * (1.0 + 2.0 * dc_cosim + di_cosisq)
                let di_f221 = 1.5 * di_sini2
                let di_f321 =  1.875 * dc_sinim * (1.0 - 2.0 * dc_cosim - 3.0 * di_cosisq)
                let di_f322 = -1.875 * dc_sinim * (1.0 + 2.0 * dc_cosim - 3.0 * di_cosisq)
                let di_f441 = 35.0 * di_sini2 * di_f220
                let di_f442 = 39.3750 * di_sini2 * di_sini2
                let di_f522 = 9.84375 * dc_sinim * (di_sini2 * (1.0 - 2.0 * dc_cosim - 5.0 * di_cosisq)
                    + 0.33333333 * (-2.0 + 4.0 * dc_cosim + 6.0 * di_cosisq))
                let di_f523 = dc_sinim * (4.92187512 * di_sini2 * (-2.0 - 4.0 * dc_cosim
                    + 10.0 * di_cosisq) + 6.56250012 * (1.0 + 2.0 * dc_cosim - 3.0 * di_cosisq))
                let di_f542 = 29.53125 * dc_sinim * (2.0 - 8.0 * dc_cosim + di_cosisq *
                    (-12.0 + 8.0 * dc_cosim + 10.0 * di_cosisq))
                let di_f543 = 29.53125 * dc_sinim * (-2.0 - 8.0 * dc_cosim + di_cosisq *
                    (12.0 + 8.0 * dc_cosim - 10.0 * di_cosisq))

                let root22: Double = 1.7891679e-6
                let root32: Double = 3.7393792e-7
                let root44: Double = 7.3636953e-9
                let root52: Double = 1.1428639e-7
                let root54: Double = 2.1765803e-9

                let di_xno2  = dc_nm * dc_nm
                let di_aonv  = pow(dc_nm / C.xke, C.x2o3)
                let di_ainv2 = di_aonv * di_aonv
                var di_temp1 = 3.0 * di_xno2 * di_ainv2
                var di_temp  = di_temp1 * root22
                ds_d2201 = di_temp * di_f220 * di_g201
                ds_d2211 = di_temp * di_f221 * di_g211
                di_temp1 = di_temp1 * di_aonv
                di_temp  = di_temp1 * root32
                ds_d3210 = di_temp * di_f321 * di_g310
                ds_d3222 = di_temp * di_f322 * di_g322
                di_temp1 = di_temp1 * di_aonv
                di_temp  = 2.0 * di_temp1 * root44
                ds_d4410 = di_temp * di_f441 * di_g410
                ds_d4422 = di_temp * di_f442 * di_g422
                di_temp1 = di_temp1 * di_aonv
                di_temp  = di_temp1 * root52
                ds_d5220 = di_temp * di_f522 * di_g520
                ds_d5232 = di_temp * di_f523 * di_g532
                di_temp  = 2.0 * di_temp1 * root54
                ds_d5421 = di_temp * di_f542 * di_g521
                ds_d5433 = di_temp * di_f543 * di_g533
                ds_xlamo = (mo + nodeo + nodeo - ds_gsto - ds_gsto)
                    .truncatingRemainder(dividingBy: C.twoPi)
                ds_xfact = mdot + ds_dmdt + 2.0 * (nodedot + ds_dnodt - rptim) - no_unkozai
                dc_em   = di_emo
                // dc_emsq restored below via dc_em
            }

            // ── Synchronous resonance (1-day orbit) ─────────────

            if ds_irez == 1 {
                let di_g200 = 1.0 + dc_emsq * (-2.5 + 0.8125 * dc_emsq)
                let di_g310 = 1.0 + 2.0 * dc_emsq
                let di_g300 = 1.0 + dc_emsq * (-6.0 + 6.60937 * dc_emsq)
                let di_f220 = 0.75 * (1.0 + dc_cosim) * (1.0 + dc_cosim)
                let di_f311 = 0.9375 * dc_sinim * dc_sinim * (1.0 + 3.0 * dc_cosim)
                    - 0.75 * (1.0 + dc_cosim)
                var di_f330 = 1.0 + dc_cosim
                di_f330 = 1.875 * di_f330 * di_f330 * di_f330

                let q22: Double = 1.7891679e-6
                let q31: Double = 2.1460748e-6
                let q33: Double = 2.2123015e-7

                let di_aonv = pow(dc_nm / C.xke, C.x2o3)
                ds_del1 = 3.0 * dc_nm * dc_nm * di_aonv * di_aonv
                ds_del2 = 2.0 * ds_del1 * di_f220 * di_g200 * q22
                ds_del3 = 3.0 * ds_del1 * di_f330 * di_g300 * q33 * di_aonv
                ds_del1 = ds_del1 * di_f311 * di_g310 * q31 * di_aonv
                ds_xlamo = (mo + nodeo + argpo - ds_gsto)
                    .truncatingRemainder(dividingBy: C.twoPi)
                ds_xfact = mdot + xpidot - rptim + ds_dmdt + ds_domdt + ds_dnodt - no_unkozai
            }

            // Initialize integrator
            ds_xli   = ds_xlamo
            ds_xni   = no_unkozai
            ds_atime = 0
            dc_nm    = no_unkozai + 0  // dndt=0 at epoch
        }

        // ── Higher-order drag (D2–D4, T3–T5) for high perigee ────

        if !isimp {
            let cc1sq = cc1 * cc1
            d2 = 4.0 * ao * tsi * cc1sq
            let temp_ = d2 * tsi * cc1 / 3.0
            d3 = (17.0 * ao + sfour) * temp_
            d4 = 0.5 * temp_ * ao * tsi * (221.0 * ao + 31.0 * sfour) * cc1
            t3cof = d2 + 2.0 * cc1sq
            t4cof = 0.25 * (3.0 * d3 + cc1 * (12.0 * d2 + 10.0 * cc1sq))
            t5cof = 0.2 * (3.0 * d4 + 12.0 * cc1 * d3
                           + 6.0 * d2 * d2 + 15.0 * cc1sq * (2.0 * d2 + cc1sq))
        } else {
            d2 = 0;  d3 = 0;  d4 = 0
            t3cof = 0;  t4cof = 0;  t5cof = 0
        }
    }

    // MARK: - Propagation

    /// Propagates to the given date and returns the ECI state vector (TEME frame).
    public nonisolated func propagate(date: Date) -> SGP4ECIState? {
        let jd = Self.julianDate(from: date)
        let tsince = (jd - epochJd) * C.minPerDay           // minutes since epoch

        // ── Secular gravity + atmospheric drag ──────────────────

        let xmdf   = mo    + mdot    * tsince
        let argpdf = argpo + argpdot * tsince
        let nodedf = nodeo + nodedot * tsince

        var argpm = argpdf
        var mm    = xmdf
        let t2    = tsince * tsince
        var nodem = nodedf + nodecf * t2
        var tempa = 1.0 - cc1 * tsince
        var tempe = bstar * cc4 * tsince
        var templ = t2cof * t2

        // ── Extended drag terms (high perigee) ──────────────────

        if !isimp {
            let delomg = omgcof * tsince
            let delmtemp = 1.0 + eta * cos(xmdf)
            let delm = xmcof * (delmtemp * delmtemp * delmtemp - delmo)
            let temp_ = delomg + delm
            mm    = xmdf + temp_
            argpm = argpdf - temp_
            let t3 = t2 * tsince
            let t4 = t3 * tsince
            tempa = tempa - d2 * t2 - d3 * t3 - d4 * t4
            tempe = tempe + bstar * cc5 * (sin(mm) - sinmao)
            templ = templ + t3cof * t3 + t4 * (t4cof + tsince * t5cof)
        }

        var nm    = no_unkozai
        var em    = ecco
        var inclm = inclo

        // Deep-space secular perturbations and resonance integration
        if isDeepSpace {
            dspace(tsince: tsince, em: &em, argpm: &argpm, inclm: &inclm,
                   mm: &mm, nodem: &nodem, nm: &nm)
        }

        // Mean motion and semi-major axis with drag correction
        let am = pow(C.xke / nm, C.x2o3) * tempa * tempa
        nm = C.xke / pow(am, 1.5)
        em = em - tempe

        // Eccentricity guard
        guard em > -0.001 else { return nil }
        if em < 1.0e-6 { em = 1.0e-6 }
        guard em < 1.0 else { return nil }

        // Mean anomaly + mean longitude
        mm = mm + no_unkozai * templ
        var xlm = mm + argpm + nodem
        nodem = nodem.truncatingRemainder(dividingBy: C.twoPi)
        argpm = argpm.truncatingRemainder(dividingBy: C.twoPi)
        xlm   = xlm.truncatingRemainder(dividingBy: C.twoPi)
        mm    = (xlm - argpm - nodem).truncatingRemainder(dividingBy: C.twoPi)

        // ── Long-period periodics ───────────────────────────────

        var ep    = em
        var xincp = inclm
        var argpp = argpm
        var nodep = nodem
        var mp    = mm

        // Deep-space periodic perturbations (dpper)
        if isDeepSpace {
            dpper(tsince: tsince, ep: &ep, inclp: &xincp,
                  nodep: &nodep, argpp: &argpp, mp: &mp)
            inclm = xincp
            if xincp < 0 {
                xincp = -xincp
                nodep += Double.pi
                argpp -= Double.pi
            }
            guard ep > 0 && ep < 1.0 else { return nil }
        }

        let sinip = sin(xincp)
        let cosip = cos(xincp)

        // For deep-space, recompute aycof/xlcof from perturbed inclination
        let aycof_: Double
        let xlcof_: Double
        if isDeepSpace {
            aycof_ = -0.5 * C.j3oj2 * sinip
            xlcof_ = abs(cosip + 1.0) > 1.5e-12
                ? -0.25 * C.j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip)
                : -0.25 * C.j3oj2 * sinip * (3.0 + 5.0 * cosip) / 1.5e-12
        } else {
            aycof_ = aycof
            xlcof_ = xlcof
        }

        let axnl = ep * cos(argpp)
        let temp_ = 1.0 / (am * (1.0 - ep * ep))
        let aynl = ep * sin(argpp) + temp_ * aycof_
        let xl   = mp + argpp + nodep + temp_ * xlcof_ * axnl

        // ── Kepler's equation ───────────────────────────────────

        var u = (xl - nodep).truncatingRemainder(dividingBy: C.twoPi)
        if u < 0 { u += C.twoPi }

        var eo1 = u
        var sineo1 = 0.0
        var coseo1 = 0.0

        for _ in 0..<10 {
            sineo1 = sin(eo1)
            coseo1 = cos(eo1)
            let denom = 1.0 - coseo1 * axnl - sineo1 * aynl
            guard abs(denom) > 1e-20 else { break }
            var tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / denom
            if abs(tem5) >= 0.95 { tem5 = tem5 > 0 ? 0.95 : -0.95 }
            eo1 += tem5
            if abs(tem5) < 1.0e-12 { break }
        }

        // ── Short-period preliminary quantities ─────────────────

        let ecose = axnl * coseo1 + aynl * sineo1
        let esine = axnl * sineo1 - aynl * coseo1
        let el2   = axnl * axnl + aynl * aynl
        let pl    = am * (1.0 - el2)
        guard pl >= 0 else { return nil }

        let rl     = am * (1.0 - ecose)
        let rdotl  = sqrt(am) * esine / rl
        let rvdotl = sqrt(pl) / rl
        let betal  = sqrt(1.0 - el2)
        let tmp    = esine / (1.0 + betal)

        let sinu = am / rl * (sineo1 - aynl - axnl * tmp)
        let cosu = am / rl * (coseo1 - axnl + aynl * tmp)
        let su   = atan2(sinu, cosu)
        let sin2u = 2.0 * cosu * sinu
        let cos2u = 1.0 - 2.0 * sinu * sinu

        let invPl  = 1.0 / pl
        let temp1  = 0.5 * C.j2 * invPl
        let temp2  = temp1 * invPl

        // ── Short-period periodics ──────────────────────────────

        let mrt = rl * (1.0 - 1.5 * temp2 * betal * con41)
                + 0.5 * temp1 * x1mth2 * cos2u
        guard mrt >= 1.0 else { return nil }   // satellite has decayed

        let suCorr   = su - 0.25 * temp2 * x7thm1 * sin2u
        let xnode    = nodep + 1.5 * temp2 * cosip * sin2u
        let xinc     = xincp + 1.5 * temp2 * cosip * sinip * cos2u
        let mvt      = rdotl  - nm * temp1 * x1mth2 * sin2u / C.xke
        let rvdot    = rvdotl + nm * temp1 * (x1mth2 * cos2u + 1.5 * con41) / C.xke

        // ── Orientation vectors → ECI (TEME) ────────────────────

        let sinsu = sin(suCorr)
        let cossu = cos(suCorr)
        let snod  = sin(xnode)
        let cnod  = cos(xnode)
        let sini  = sin(xinc)
        let cosi  = cos(xinc)

        let xmx = -snod * cosi
        let xmy =  cnod * cosi

        let ux = xmx * sinsu + cnod * cossu
        let uy = xmy * sinsu + snod * cossu
        let uz = sini * sinsu

        let vx = xmx * cossu - cnod * sinsu
        let vy = xmy * cossu - snod * sinsu
        let vz = sini * cossu

        // Position (km) and velocity (km/s)
        let position = SIMD3<Double>(
            mrt * ux * C.earthRadius,
            mrt * uy * C.earthRadius,
            mrt * uz * C.earthRadius
        )
        let velocity = SIMD3<Double>(
            (mvt * ux + rvdot * vx) * C.vkmpersec,
            (mvt * uy + rvdot * vy) * C.vkmpersec,
            (mvt * uz + rvdot * vz) * C.vkmpersec
        )

        return SGP4ECIState(position: position, velocity: velocity)
    }

    // MARK: - Deep-Space Perturbation Methods

    /// Deep-space periodic perturbations (solar/lunar).
    /// Called during propagation for deep-space orbits (period ≥ 225 min).
    nonisolated private func dpper(
        tsince: Double,
        ep: inout Double, inclp: inout Double,
        nodep: inout Double, argpp: inout Double, mp: inout Double
    ) {
        let zes: Double = 0.01675
        let zel: Double = 0.05490
        let zns: Double = 1.19459e-5
        let znl: Double = 1.5835218e-4

        // Solar terms
        var zm = ds_zmos + zns * tsince
        var zf = zm + 2.0 * zes * sin(zm)
        var sinzf = sin(zf)
        var f2 = 0.5 * sinzf * sinzf - 0.25
        var f3 = -0.5 * sinzf * cos(zf)
        let ses  = ds_se2 * f2 + ds_se3 * f3
        let sis  = ds_si2 * f2 + ds_si3 * f3
        let sls  = ds_sl2 * f2 + ds_sl3 * f3 + ds_sl4 * sinzf
        let sghs = ds_sgh2 * f2 + ds_sgh3 * f3 + ds_sgh4 * sinzf
        let shs  = ds_sh2 * f2 + ds_sh3 * f3

        // Lunar terms
        zm = ds_zmol + znl * tsince
        zf = zm + 2.0 * zel * sin(zm)
        sinzf = sin(zf)
        f2 = 0.5 * sinzf * sinzf - 0.25
        f3 = -0.5 * sinzf * cos(zf)
        let sel  = ds_ee2 * f2 + ds_e3 * f3
        let sil  = ds_xi2 * f2 + ds_xi3 * f3
        let sll  = ds_xl2 * f2 + ds_xl3 * f3 + ds_xl4 * sinzf
        let sghl = ds_xgh2 * f2 + ds_xgh3 * f3 + ds_xgh4 * sinzf
        let shll = ds_xh2 * f2 + ds_xh3 * f3

        // Sum solar + lunar, subtract epoch values (relative correction)
        let pe   = ses + sel - ds_peo
        let pinc = sis + sil - ds_pinco
        let pl   = sls + sll - ds_plo
        var pgh  = sghs + sghl - ds_pgho
        var ph   = shs + shll - ds_pho

        inclp += pinc
        ep    += pe
        let sinip = sin(inclp)
        let cosip = cos(inclp)

        if inclp >= 0.2 {
            ph  = ph / sinip
            pgh = pgh - cosip * ph
            argpp += pgh
            nodep += ph
            mp    += pl
        } else {
            let sinop = sin(nodep)
            let cosop = cos(nodep)
            var alfdp = sinip * sinop
            var betdp = sinip * cosop
            let dalf  =  ph * cosop + pinc * cosip * sinop
            let dbet  = -ph * sinop + pinc * cosip * cosop
            alfdp += dalf
            betdp += dbet
            nodep = nodep.truncatingRemainder(dividingBy: C.twoPi)
            if nodep < 0 { nodep += C.twoPi }
            let xls = mp + argpp + pl + pgh
                + (cosip - pinc * sinip) * nodep
            let xnoh = nodep
            nodep = atan2(alfdp, betdp)
            if abs(xnoh - nodep) > Double.pi {
                nodep = nodep < xnoh ? nodep + C.twoPi : nodep - C.twoPi
            }
            mp += pl
            argpp = xls - mp - cosip * nodep
        }
    }

    /// Deep-space secular effects and resonance integration.
    /// Updates mean elements for solar/lunar gravity and performs
    /// Euler-Maclaurin numerical integration for resonant orbits.
    nonisolated private func dspace(
        tsince: Double,
        em: inout Double, argpm: inout Double, inclm: inout Double,
        mm: inout Double, nodem: inout Double, nm: inout Double
    ) {
        let rptim: Double = 4.37526908801129966e-3
        let fasx2: Double = 0.13130908
        let fasx4: Double = 2.8843198
        let fasx6: Double = 0.37448087
        let g22: Double   = 5.7686396
        let g32: Double   = 0.95240898
        let g44: Double   = 1.8014998
        let g52: Double   = 1.0508330
        let g54: Double   = 4.4108898
        let stepp: Double =    720.0
        let stepn: Double =   -720.0
        let step2: Double = 259200.0

        let theta = (ds_gsto + tsince * rptim)
            .truncatingRemainder(dividingBy: C.twoPi)

        // Secular rates from lunar/solar gravity
        em    += ds_dedt * tsince
        inclm += ds_didt * tsince
        argpm += ds_domdt * tsince
        nodem += ds_dnodt * tsince
        mm    += ds_dmdt * tsince

        // Resonance integration (Euler-Maclaurin)
        if ds_irez != 0 {
            // Epoch restart if needed
            if ds_atime == 0 || tsince * ds_atime <= 0 || abs(tsince) < abs(ds_atime) {
                ds_atime = 0
                ds_xni = no_unkozai
                ds_xli = ds_xlamo
            }

            let delt: Double = tsince > 0 ? stepp : stepn
            var ft = 0.0
            var xndt = 0.0, xldot = 0.0, xnddt = 0.0

            var stepping = true
            while stepping {
                // Compute derivatives at current integration point
                if ds_irez != 2 {
                    // Synchronous resonance (1-day orbit)
                    xndt = ds_del1 * sin(ds_xli - fasx2)
                        + ds_del2 * sin(2.0 * (ds_xli - fasx4))
                        + ds_del3 * sin(3.0 * (ds_xli - fasx6))
                    xldot = ds_xni + ds_xfact
                    xnddt = ds_del1 * cos(ds_xli - fasx2)
                        + 2.0 * ds_del2 * cos(2.0 * (ds_xli - fasx4))
                        + 3.0 * ds_del3 * cos(3.0 * (ds_xli - fasx6))
                    xnddt *= xldot
                } else {
                    // Half-day resonance (12-hour orbit)
                    let xomi  = argpo + argpdot * ds_atime
                    let x2omi = xomi + xomi
                    let x2li  = ds_xli + ds_xli
                    xndt = ds_d2201 * sin(x2omi + ds_xli - g22)
                        + ds_d2211 * sin(ds_xli - g22)
                        + ds_d3210 * sin(xomi + ds_xli - g32)
                        + ds_d3222 * sin(-xomi + ds_xli - g32)
                        + ds_d4410 * sin(x2omi + x2li - g44)
                        + ds_d4422 * sin(x2li - g44)
                        + ds_d5220 * sin(xomi + ds_xli - g52)
                        + ds_d5232 * sin(-xomi + ds_xli - g52)
                        + ds_d5421 * sin(xomi + x2li - g54)
                        + ds_d5433 * sin(-xomi + x2li - g54)
                    xldot = ds_xni + ds_xfact
                    xnddt = ds_d2201 * cos(x2omi + ds_xli - g22)
                        + ds_d2211 * cos(ds_xli - g22)
                        + ds_d3210 * cos(xomi + ds_xli - g32)
                        + ds_d3222 * cos(-xomi + ds_xli - g32)
                        + ds_d5220 * cos(xomi + ds_xli - g52)
                        + ds_d5232 * cos(-xomi + ds_xli - g52)
                    xnddt += 2.0 * (
                        ds_d4410 * cos(x2omi + x2li - g44)
                        + ds_d4422 * cos(x2li - g44)
                        + ds_d5421 * cos(xomi + x2li - g54)
                        + ds_d5433 * cos(-xomi + x2li - g54))
                    xnddt *= xldot
                }

                // Check if we need to keep stepping
                if abs(tsince - ds_atime) >= stepp {
                    // Full integration step
                    ds_xli  += xldot * delt + xndt * step2
                    ds_xni  += xndt * delt + xnddt * step2
                    ds_atime += delt
                } else {
                    // Final fractional step
                    ft = tsince - ds_atime
                    stepping = false
                }
            }

            // Interpolate to exact time
            nm = ds_xni + xndt * ft + xnddt * ft * ft * 0.5
            let xl = ds_xli + xldot * ft + xndt * ft * ft * 0.5

            if ds_irez != 1 {
                mm = xl - 2.0 * nodem + 2.0 * theta
            } else {
                mm = xl - nodem - argpm + theta
            }
            nm = no_unkozai + (nm - no_unkozai)
        }
    }

    // MARK: - ECI → Geodetic

    /// Converts ECI position (TEME) to geodetic lat/lon/alt.
    /// Uses iterative geodetic latitude for WGS-72 oblate spheroid.
    public nonisolated static func eciToGeodetic(position: SIMD3<Double>, date: Date) -> SGP4GeodeticPosition {
        let jd = julianDate(from: date)
        let theta = gmst(jd: jd)

        let x = position.x
        let y = position.y
        let z = position.z

        // Longitude: rotate from TEME to PEF
        var lon = atan2(y, x) - theta
        // Normalize to −π … +π
        lon = lon.truncatingRemainder(dividingBy: C.twoPi)
        if lon < -Double.pi { lon += C.twoPi }
        if lon >  Double.pi { lon -= C.twoPi }

        let rxy = sqrt(x * x + y * y)

        // Iterative geodetic latitude (WGS-72 ellipsoid)
        let f  = 1.0 / 298.26          // WGS-72 flattening
        let e2 = 2.0 * f - f * f       // first eccentricity squared
        var lat = atan2(z, rxy)         // initial geocentric guess

        for _ in 0..<10 {
            let sinLat = sin(lat)
            let cRad = 1.0 / sqrt(1.0 - e2 * sinLat * sinLat)
            lat = atan2(z + C.earthRadius * e2 * cRad * sinLat, rxy)
        }

        // Altitude
        let sinLat = sin(lat)
        let cRad = 1.0 / sqrt(1.0 - e2 * sinLat * sinLat)
        let alt: Double
        if abs(lat) < Double.pi / 4.0 {
            // Near equator: use rxy
            alt = rxy / cos(lat) - C.earthRadius * cRad
        } else {
            // Near poles: use z
            alt = z / sinLat - C.earthRadius * cRad * (1.0 - e2)
        }

        return SGP4GeodeticPosition(
            latitude:  lat * C.rad2deg,
            longitude: lon * C.rad2deg,
            altitude:  alt
        )
    }

    // MARK: - Look Angles

    /// Computes azimuth, elevation, range, and range-rate from an observer to the satellite.
    public nonisolated static func lookAngles(
        position: SIMD3<Double>,
        velocity: SIMD3<Double>,
        observer: CLLocation,
        date: Date
    ) -> SGP4LookAngles {
        let jd = julianDate(from: date)
        let theta = gmst(jd: jd)

        let latRad = observer.coordinate.latitude  * C.deg2rad
        let lonRad = observer.coordinate.longitude * C.deg2rad
        let altKm  = observer.altitude / 1000.0

        let obsTheta = theta + lonRad

        // Observer position in ECI (spherical approximation, adequate for tracking)
        let cosLat   = cos(latRad)
        let sinLat   = sin(latRad)
        let cosTheta = cos(obsTheta)
        let sinTheta = sin(obsTheta)

        let obsR  = C.earthRadius + altKm
        let obsX  = obsR * cosLat * cosTheta
        let obsY  = obsR * cosLat * sinTheta
        let obsZ  = obsR * sinLat

        // Observer velocity in ECI (Earth rotation)
        let obsVx = -C.earthRotRate * obsY
        let obsVy =  C.earthRotRate * obsX
        let obsVz =  0.0

        // Range vector (ECI)
        let rx = position.x - obsX
        let ry = position.y - obsY
        let rz = position.z - obsZ
        let range = sqrt(rx * rx + ry * ry + rz * rz)

        // Relative velocity
        let drx = velocity.x - obsVx
        let dry = velocity.y - obsVy
        let drz = velocity.z - obsVz
        let rangeRate = (rx * drx + ry * dry + rz * drz) / max(range, 1e-10)

        // Topocentric (SEZ) rotation
        let topS =  sinLat * cosTheta * rx + sinLat * sinTheta * ry - cosLat * rz
        let topE = -sinTheta * rx + cosTheta * ry
        let topZ =  cosLat * cosTheta * rx + cosLat * sinTheta * ry + sinLat * rz

        var azimuth = atan2(-topE, topS) + Double.pi
        azimuth = azimuth.truncatingRemainder(dividingBy: C.twoPi)
        if azimuth < 0 { azimuth += C.twoPi }

        let elevation = asin(min(1.0, max(-1.0, topZ / max(range, 1e-10))))

        return SGP4LookAngles(
            azimuth:   azimuth   * C.rad2deg,
            elevation: elevation * C.rad2deg,
            range:     range,
            rangeRate: rangeRate
        )
    }

    // MARK: - Time Utilities

    /// Converts a Swift Date to Julian Date.
    public nonisolated static func julianDate(from date: Date) -> Double {
        let cal = Calendar(identifier: .gregorian)
        let c   = cal.dateComponents(in: TimeZone(identifier: "UTC")!, from: date)
        let yr  = Double(c.year!)
        let mon = Double(c.month!)
        let day = Double(c.day!)
        let hr  = Double(c.hour   ?? 0)
        let min = Double(c.minute ?? 0)
        let sec = Double(c.second ?? 0)

        var y = yr, m = mon
        if mon <= 2 { y -= 1; m += 12 }

        let A = floor(y / 100.0)
        let B = 2.0 - A + floor(A / 4.0)
        let jdDay = floor(365.25 * (y + 4716.0)) + floor(30.6001 * (m + 1.0)) + day + B - 1524.5
        let jdFrac = (hr + min / 60.0 + sec / 3600.0) / 24.0

        return jdDay + jdFrac
    }

    /// Greenwich Mean Sidereal Time in radians.
    /// Uses the IAU 1982 / Vallado formulation.
    public nonisolated static func gmst(jd: Double) -> Double {
        let tut1 = (jd - 2451545.0) / 36525.0
        var temp = -6.2e-6 * tut1 * tut1 * tut1
                 + 0.093104 * tut1 * tut1
                 + (876600.0 * 3600.0 + 8640184.812866) * tut1
                 + 67310.54841                               // seconds
        temp = (temp * C.deg2rad / 240.0)                      // → radians
            .truncatingRemainder(dividingBy: C.twoPi)
        if temp < 0 { temp += C.twoPi }
        return temp
    }

    // MARK: - Private Helpers

    /// Julian Date from TLE epoch year + fractional day of year.
    nonisolated private static func julianDateFromEpoch(year: Int, day: Double) -> Double {
        let yr = Double(year)
        let base = 367.0 * yr
                 - floor(7.0 * (yr + floor(10.0 / 12.0)) / 4.0)
                 + floor(275.0 / 9.0)
                 + 1721013.5
        return base + day
    }

    /// Parses TLE exponential notation (e.g., "12345-4" → 0.12345×10⁻⁴).
    nonisolated private static func parseExponentialField(_ field: String) -> Double {
        var s = field.trimmingCharacters(in: .whitespaces)
        guard !s.isEmpty else { return 0 }

        var sign = 1.0
        if s.first == "-" { sign = -1.0; s = String(s.dropFirst()) }
        else if s.first == "+" { s = String(s.dropFirst()) }

        if let expIdx = s.lastIndex(where: { $0 == "-" || $0 == "+" }) {
            let mantissa = String(s[s.startIndex..<expIdx])
            let exponent = String(s[expIdx...])
            if let m = Double("0." + mantissa), let e = Double(exponent) {
                return sign * m * pow(10.0, e)
            }
        }

        return Double(s).map { sign * $0 } ?? 0
    }
}
