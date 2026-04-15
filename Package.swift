// swift-tools-version: 5.9

import PackageDescription

let package = Package(
    name: "AllMySatKit",
    platforms: [
        .iOS(.v14),
        .macOS(.v12),
        .watchOS(.v7),
        .tvOS(.v14),
        .visionOS(.v1)
    ],
    products: [
        .library(name: "AllMySatKit", targets: ["AllMySatKit"]),
    ],
    targets: [
        .target(
            name: "AllMySatKit",
            path: "Sources/AllMySatKit",
            exclude: [
                "Propagation/README.md",
                "Astronomy/README.md",
                "Geodesy/README.md",
                "TLE/README.md",
                "Radio/README.md"
            ]
        ),
        .testTarget(
            name: "AllMySatKitTests",
            dependencies: ["AllMySatKit"],
            path: "Tests/AllMySatKitTests"
        ),
    ]
)
