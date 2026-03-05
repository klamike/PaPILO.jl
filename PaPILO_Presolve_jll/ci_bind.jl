#!/usr/bin/env julia
# Called from the release job of build_jll.yml.
# Reads every artifact_meta_*.txt in the current directory, then updates
# Artifacts.toml to point to the freshly created GitHub release.
#
# Env vars required:
#   REPO           – "owner/repo", e.g. "klamike/PaPILO.jl"
#   TAG            – release tag, e.g. "PaPILO_Presolve-v3.1.0"
#   ARTIFACT_NAME  – name key in Artifacts.toml, e.g. "PaPILO_Presolve"
#   ARTIFACTS_TOML – absolute path to the Artifacts.toml to update

using Pkg.Artifacts
using Base.BinaryPlatforms

repo      = ENV["REPO"]
tag       = ENV["TAG"]
art_name  = ENV["ARTIFACT_NAME"]
toml_path = ENV["ARTIFACTS_TOML"]

meta_files = filter(readdir(".")) do f
    startswith(f, "artifact_meta_") && endswith(f, ".txt")
end
isempty(meta_files) && error("No artifact_meta_*.txt files found in $(pwd())")

for f in meta_files
    d = Dict{String,String}()
    for l in readlines(f)
        contains(l, "=") || continue
        k, v = split(l, "=", limit=2)
        d[k] = v
    end

    platform_str = d["PLATFORM"]
    tree_sha1    = d["TREE_SHA1"]
    sha256_hex   = d["SHA256"]
    tarball_name = d["TARBALL"]

    url = "https://github.com/$(repo)/releases/download/$(tag)/$(tarball_name)"

    arch = String(split(platform_str, "-")[1])
    plat = if contains(platform_str, "linux")
        Platform(arch, "linux")      # libc=glibc by default
    elseif contains(platform_str, "apple")
        Platform(arch, "macos")
    else
        error("Unrecognised platform string: $platform_str")
    end

    tree_hash = Base.SHA1(hex2bytes(tree_sha1))

    bind_artifact!(toml_path, art_name, tree_hash;
        platform      = plat,
        download_info = [(url, sha256_hex)],
        force         = true,
    )
    @info "Bound $(platform_str)" url
end
