#!/usr/bin/env julia
# Build libpapilo-presolve.dylib locally and repack the JLL artifact.
# Requires an already-configured cmake build at PAPILO_BUILD_DIR.
# Run from anywhere:  julia /path/to/build_local.jl

using Pkg.Artifacts
using SHA: sha256

const JLLDIR       = @__DIR__
const REPO         = normpath(joinpath(JLLDIR, "../.."))
const PAPILO_SRC   = joinpath(REPO, "papilo")
const PAPILO_BUILD = joinpath(PAPILO_SRC, "build")
const BUNDLED_SRC  = joinpath(JLLDIR, "bundled", "src")
const PATCH        = joinpath(JLLDIR, "bundled", "patches", "papilo_presolve.patch")
const PRODUCTS     = joinpath(JLLDIR, "products")
const ARTIFACTS_FILE = joinpath(JLLDIR, "Artifacts.toml")

mkpath(PRODUCTS)

isdir(PAPILO_BUILD) || error(
    "PaPILO build directory not found: $PAPILO_BUILD\n" *
    "Configure PaPILO first with cmake before running this script.")

# ── 0. Inject bundled sources and apply patch ──────────────────────────────────
@info "Copying bundled wrapper sources to papilo/src/..."
for f in ("papilo_presolve.h", "papilo_presolve.cpp")
    cp(joinpath(BUNDLED_SRC, f), joinpath(PAPILO_SRC, "src", f); force = true)
end

@info "Applying CMakeLists patch..."
run(`git -C $PAPILO_SRC apply $PATCH`)

# ── 1. Re-run cmake to pick up the new papilo-presolve target, then build ──────
@info "Reconfiguring cmake..."
tbb_dir = "/opt/homebrew/opt/tbb/lib/cmake/TBB"
run(`cmake -S $PAPILO_SRC -B $PAPILO_BUILD
    -DCMAKE_BUILD_TYPE=Release
    -DBUILD_TESTING=OFF
    -DPAPILO_NO_BINARIES=ON
    -DTBB=ON
    -DTBB_DIR=$tbb_dir`)

@info "Building libpapilo-presolve..."
run(`cmake --build $PAPILO_BUILD --target papilo-presolve -j$(Sys.CPU_THREADS)`)

# ── 2. Locate and copy the built library ──────────────────────────────────────
libname  = Sys.isapple() ? "libpapilo-presolve.dylib" : "libpapilo-presolve.so"
built_lib = joinpath(PAPILO_BUILD, libname)
isfile(built_lib) || error("Expected $built_lib not found after build")
out_lib = joinpath(PRODUCTS, libname)
cp(built_lib, out_lib; force = true)
@info "Built: $out_lib"

# ── 3. Create the artifact ─────────────────────────────────────────────────────
@info "Creating artifact..."
tree_hash = create_artifact() do dir
    libdir = joinpath(dir, "lib")
    mkpath(libdir)
    cp(out_lib, joinpath(libdir, libname))
end

# ── 4. Archive to products/ ────────────────────────────────────────────────────
platform = "aarch64-apple-darwin-libgfortran5"
version  = "3.1.0"
tarball  = joinpath(PRODUCTS, "PaPILO_Presolve.v$(version).$(platform).tar.gz")
archive_artifact(tree_hash, tarball)

sha256_hex = bytes2hex(open(sha256, tarball))
url = "file://$(tarball)"
@info "Tarball: $tarball (sha256=$sha256_hex)"

# ── 5. Bind into Artifacts.toml ────────────────────────────────────────────────
plat = Base.BinaryPlatforms.Platform("aarch64", "macos";
    libgfortran_version = "5.0.0")
bind_artifact!(ARTIFACTS_FILE, "PaPILO_Presolve", tree_hash;
    platform      = plat,
    download_info = [(url, sha256_hex)],
    force         = true,
)
@info "Artifacts.toml updated."

# ── 6. Revert patch (keep papilo source clean) ─────────────────────────────────
@info "Reverting patch..."
run(`git -C $PAPILO_SRC apply -R $PATCH`)
# Remove injected files
for f in ("papilo_presolve.h", "papilo_presolve.cpp")
    rm(joinpath(PAPILO_SRC, "src", f); force = true)
end
@info "Done."
