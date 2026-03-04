module API

using PaPILO_Presolve_jll: libpapilo_presolve

@enum Status begin
    Unchanged             = 0
    Reduced               = 1
    Infeasible            = 2
    Unbounded             = 3
    UnboundedOrInfeasible = 4
end

mutable struct Presolver
    _ptr::Ptr{Cvoid}

    function Presolver(; threads::Int = 1, verbosity::Int = 0)
        ptr = ccall((:papilo_create, libpapilo_presolve), Ptr{Cvoid}, ())
        ptr == C_NULL && error("papilo_create failed: $(_last_error())")
        obj = new(ptr)
        finalizer(p -> ccall((:papilo_free, libpapilo_presolve), Cvoid, (Ptr{Cvoid},), p._ptr), obj)
        ccall((:papilo_set_threads,   libpapilo_presolve), Cvoid, (Ptr{Cvoid}, Cint), ptr, threads)
        ccall((:papilo_set_verbosity, libpapilo_presolve), Cvoid, (Ptr{Cvoid}, Cint), ptr, verbosity)
        return obj
    end
end

mutable struct Result
    _ptr::Ptr{Cvoid}

    function Result(ptr::Ptr{Cvoid})
        obj = new(ptr)
        finalizer(r -> ccall((:papilo_result_free, libpapilo_presolve), Cvoid, (Ptr{Cvoid},), r._ptr), obj)
        return obj
    end
end

_last_error() = begin
    ptr = ccall((:papilo_last_error, libpapilo_presolve), Ptr{UInt8}, ())
    ptr == C_NULL ? "(no error)" : unsafe_string(ptr)
end

_inf_flags(v::AbstractVector{Float64}) = UInt8.(isinf.(v))
_strip_inf(v::AbstractVector{Float64}) = map(x -> isinf(x) ? 0.0 : x, v)

function _restore_inf(v::Vector{Float64}, flags::Vector{UInt8}, negative::Bool)
    fill_val = negative ? -Inf : Inf
    return [flags[i] != 0 ? fill_val : v[i] for i in eachindex(v)]
end

function apply(
        p          :: Presolver,
        ncols      :: Int,
        nrows      :: Int,
        c          :: Vector{Float64},
        obj_offset :: Float64,
        col_lb     :: Vector{Float64},
        col_ub     :: Vector{Float64},
        row_lhs    :: Vector{Float64},
        row_rhs    :: Vector{Float64},
        row_start  :: Vector{Int},   # 1-indexed Julia colptr / rowptr
        col_idx    :: Vector{Int},   # 1-indexed Julia rowval / colval
        vals       :: Vector{Float64},
    )
    rs       = Int32.(row_start .- 1)   # 1-indexed → 0-indexed
    ci       = Int32.(col_idx  .- 1)   # 1-indexed → 0-indexed
    lbi      = _inf_flags(col_lb);  ubi  = _inf_flags(col_ub)
    lhi      = _inf_flags(row_lhs); rhi  = _inf_flags(row_rhs)
    integral = zeros(UInt8, ncols)   # treat all variables as continuous

    ptr = ccall(
        (:papilo_apply, libpapilo_presolve), Ptr{Cvoid},
        (
            Ptr{Cvoid}, Cint, Cint,
            Ptr{Cdouble}, Cdouble,
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cuchar}, Ptr{Cuchar}, Ptr{Cuchar},
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cuchar}, Ptr{Cuchar},
            Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
        ),
        p._ptr, Cint(ncols), Cint(nrows),
        _strip_inf(c), Cdouble(obj_offset),
        _strip_inf(col_lb), _strip_inf(col_ub), lbi, ubi, integral,
        _strip_inf(row_lhs), _strip_inf(row_rhs), lhi, rhi,
        Cint(length(vals)), rs, ci, vals,
    )
    ptr == C_NULL && error("papilo_apply failed: $(_last_error())")
    return Result(ptr)
end

status(r::Result)    = Status(Int(ccall((:papilo_status,   libpapilo_presolve), Cint, (Ptr{Cvoid},), r._ptr)))
num_cols(r::Result)  = Int(ccall((:papilo_num_cols,  libpapilo_presolve), Cint, (Ptr{Cvoid},), r._ptr))
num_rows(r::Result)  = Int(ccall((:papilo_num_rows,  libpapilo_presolve), Cint, (Ptr{Cvoid},), r._ptr))
nnz(r::Result)       = Int(ccall((:papilo_nnz,       libpapilo_presolve), Cint, (Ptr{Cvoid},), r._ptr))
orig_cols(r::Result) = Int(ccall((:papilo_orig_cols, libpapilo_presolve), Cint, (Ptr{Cvoid},), r._ptr))
orig_rows(r::Result) = Int(ccall((:papilo_orig_rows, libpapilo_presolve), Cint, (Ptr{Cvoid},), r._ptr))

function get_obj(r::Result)
    coeffs = Vector{Float64}(undef, num_cols(r))
    offset = Ref{Float64}(0.0)
    ccall((:papilo_get_obj, libpapilo_presolve), Cvoid,
          (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), r._ptr, coeffs, offset)
    return coeffs, offset[]
end

function get_col_bounds(r::Result)
    n   = num_cols(r)
    lb  = Vector{Float64}(undef, n);  ub  = Vector{Float64}(undef, n)
    lbi = Vector{UInt8}(undef, n);    ubi = Vector{UInt8}(undef, n)
    iti = Vector{UInt8}(undef, n)
    ccall((:papilo_get_col_bounds, libpapilo_presolve), Cvoid,
          (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cuchar}, Ptr{Cuchar}, Ptr{Cuchar}),
          r._ptr, lb, ub, lbi, ubi, iti)
    return _restore_inf(lb, lbi, true), _restore_inf(ub, ubi, false)
end

function get_row_bounds(r::Result)
    m   = num_rows(r)
    lhs = Vector{Float64}(undef, m);  rhs = Vector{Float64}(undef, m)
    lhi = Vector{UInt8}(undef, m);    rhi = Vector{UInt8}(undef, m)
    ccall((:papilo_get_row_bounds, libpapilo_presolve), Cvoid,
          (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cuchar}, Ptr{Cuchar}),
          r._ptr, lhs, rhs, lhi, rhi)
    return _restore_inf(lhs, lhi, true), _restore_inf(rhs, rhi, false)
end

function get_matrix(r::Result)
    m  = num_rows(r);  nz = nnz(r)
    rs  = Vector{Int32}(undef, m + 1)
    ci  = Vector{Int32}(undef, nz)
    nzv = Vector{Float64}(undef, nz)
    ccall((:papilo_get_matrix, libpapilo_presolve), Cvoid,
          (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
          r._ptr, rs, ci, nzv)
    return Vector{Int}(rs) .+ 1, Vector{Int}(ci) .+ 1, nzv   # 0-indexed → 1-indexed
end

function get_col_map(r::Result)
    buf = Vector{Int32}(undef, num_cols(r))
    ccall((:papilo_get_col_map, libpapilo_presolve), Cvoid, (Ptr{Cvoid}, Ptr{Cint}), r._ptr, buf)
    return Vector{Int}(buf) .+ 1   # 0-indexed → 1-indexed
end

function get_row_map(r::Result)
    buf = Vector{Int32}(undef, num_rows(r))
    ccall((:papilo_get_row_map, libpapilo_presolve), Cvoid, (Ptr{Cvoid}, Ptr{Cint}), r._ptr, buf)
    return Vector{Int}(buf) .+ 1   # 0-indexed → 1-indexed
end

function postsolve(r::Result, reduced::Vector{Float64})
    original = Vector{Float64}(undef, orig_cols(r))
    ret = ccall((:papilo_postsolve, libpapilo_presolve), Cint,
                (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}), r._ptr, reduced, original)
    ret != 0 && @warn "papilo_postsolve returned non-zero status" ret
    return original
end

function map_primal(r::Result, x_orig::Vector{Float64})
    x_red = Vector{Float64}(undef, num_cols(r))
    ccall((:papilo_map_primal, libpapilo_presolve), Cvoid,
          (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}),
          r._ptr, x_orig, x_red)
    return x_red
end

end  # module API
