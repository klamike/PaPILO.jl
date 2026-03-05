using BinaryBuilder, Pkg

include("macos_sdks.jl")

name = "PaPILO_Presolve"
version = v"3.0.0"

sources = [
    GitSource("https://github.com/scipopt/papilo.git", "4cbdc3270e5f454c0b6ce81454d7ee5acefa7777"),
    DirectorySource("./bundled"),
]

preamble = if os() === "macos"  # https://github.com/JuliaPackaging/BinaryBuilder.jl/issues/1348
    raw"find  /usr/share/cmake/Modules/Compiler/ -name '._*' -delete ; "
else
    raw""
end 

script = preamble * raw"""
cd ${WORKSPACE}/srcdir/papilo
cp ${WORKSPACE}/srcdir/src/papilo_presolve.h src/
cp ${WORKSPACE}/srcdir/src/papilo_presolve.cpp src/
atomic_patch -p1 ${WORKSPACE}/srcdir/patches/papilo_presolve.patch

cmake -S . -B build \
    -DCMAKE_INSTALL_PREFIX=${prefix} \
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=OFF \
    -DPAPILO_NO_BINARIES=ON \
    -DTBB=ON \
    -DTBB_DIR=${prefix}/lib/cmake/TBB \
    -DGMP=ON

cmake --build build --target papilo-presolve -j${nproc}
cmake --install build
install_license ${WORKSPACE}/srcdir/papilo/LICENSE
"""

platforms = supported_platforms()
platforms = expand_cxxstring_abis(platforms)
platforms = expand_gfortran_versions(platforms)

filter!(p -> !(Sys.isfreebsd(p) && arch(p) == "aarch64"), platforms)
filter!(p -> !(Sys.islinux(p) && arch(p) == "riscv64"), platforms)
filter!(p -> nbits(p) != 32, platforms)
filter!(p -> libgfortran_version(p) >= v"5", platforms)

products = [
    LibraryProduct("libpapilo-presolve", :libpapilo_presolve),
]

dependencies = [
    Dependency("boost_jll"; compat="=1.79.0"),
    Dependency("oneTBB_jll"; compat="2021.5.0"),
    Dependency("GMP_jll"; compat="6.2.1"),
    Dependency("CompilerSupportLibraries_jll"),
]

sources, script = require_macos_sdk("10.13", sources, script)
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies;
               julia_compat="1.6", preferred_gcc_version=v"12", clang_use_lld=false)
