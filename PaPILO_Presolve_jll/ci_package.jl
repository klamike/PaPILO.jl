#!/usr/bin/env julia
# Called from the build job of build_jll.yml.
# Env vars required:
#   VERSION       – JLL version string, e.g. "3.1.0"
#   PLATFORM      – BinaryPlatforms triplet, e.g. "x86_64-linux-gnu"
#   ARTIFACT_NAME – name used in Artifacts.toml, e.g. "PaPILO_Presolve"
#   BUILT_LIB     – absolute path to the compiled shared library

using Pkg.Artifacts, SHA

version   = ENV["VERSION"]
platform  = ENV["PLATFORM"]
art_name  = ENV["ARTIFACT_NAME"]
built     = ENV["BUILT_LIB"]

isfile(built) || error("Built library not found: $built")
libname = basename(built)

tree_hash = create_artifact() do dir
    mkpath(joinpath(dir, "lib"))
    cp(built, joinpath(dir, "lib", libname))
end

tarball = "$(art_name).v$(version).$(platform).tar.gz"
archive_artifact(tree_hash, tarball)

sha256_hex = bytes2hex(open(sha256, tarball))
tree_sha1  = string(tree_hash)

open("artifact_meta_$(platform).txt", "w") do f
    println(f, "TREE_SHA1=$(tree_sha1)")
    println(f, "SHA256=$(sha256_hex)")
    println(f, "TARBALL=$(tarball)")
    println(f, "PLATFORM=$(platform)")
end

@info "Packaged $(tarball)" tree=tree_sha1 sha256=sha256_hex
