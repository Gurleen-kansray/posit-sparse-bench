using UniversalNumbers
using SparseArrays
using LinearAlgebra
using MatrixDepot
using BPosits
using Random
using Statistics
using Printf

# 
using UniversalNumbers: Quire

const MAXITER = 300
const RTOL    = 1e-10          # ||r|| / ||b|| stopping tolerance

const DRIVEN_MAX_NNZ    = 150_000   # stage 2 budget
const SEED_MAX_N        = 5_000     # stage 3 ladder budget (ladder is O(n))
const SEED_DRIVEN_MAX_NNZ = 150_000 # stage 3 driven budget
const NSEEDS_S4_LARGE = 20
const S4_LARGE_NNZ    = 30_000
const NSEEDS            = 50        # matches the 50-seeds-per-matrix study
const SPECTRUM_MAX_N    = 2_000     # dense eigvals is O(n^3)
const SYMSCALE = false
const MATRICES = [
    ("bcsstk03", "HB/bcsstk03"),          #   112,     640
    ("nos2",     "HB/nos2"),              #   957,    4137
    ("mhd4800b", "Bai/mhd4800b"),         #  4800,   27520
    ("bcsstk14", "HB/bcsstk14"),          #  1806,   63454
    ("sts4098",  "Cannizzo/sts4098"),     #  4098,   72356
    ("nasa4704", "Nasa/nasa4704"),        #  4704,  104756
    ("bodyy4",   "Pothen/bodyy4"),        # 17546,  121938
    ("bcsstk38", "Boeing/bcsstk38"),      #  8032,  355460
    ("bcsstk37", "Boeing/bcsstk37"),      # 25503, 1140977
    ("bcsstk36", "Boeing/bcsstk36"),      # 23052, 1143140
    ("nasasrb",  "Nasa/nasasrb"),         # 54870, 2677324
    ("s3dkt3m2", "GHS_psdef/s3dkt3m2"),   # 90449, 3753461
    ("s3dkq4m2", "GHS_psdef/s3dkq4m2"),   # 90449, 4820891
]

const LADDER = ((Posit{8,2,UInt8},   BPosit8,  "p8"),
                (Posit{16,2,UInt16}, BPosit16, "p16"),
                (Posit{32,2,UInt32}, BPosit32, "p32"))



function naive_dot(a::AbstractVector, b::AbstractVector)
    s = zero(eltype(a))
    for i in eachindex(a, b)
        s += a[i] * b[i]
    end
    return s
end

quire_dot(a::AbstractVector{T}, b::AbstractVector{T}) where {T<:Posit} = fdp(a, b)
quire_dot(a::AbstractVector{Float64}, b::AbstractVector{Float64}) = naive_dot(a, b)
quire_dot(a::AbstractVector{T}, b::AbstractVector{T}) where {T<:AnyBPosit} =
    LinearAlgebra.dot(a, b; quire=true)

relerr(v, ref) = iszero(ref) ? (iszero(v) ? 0.0 : Inf) : abs(v - ref) / abs(ref)
same(x, y) = (isnan(x) && isnan(y)) || x == y

plain_matvec(A, x) = A * x

function quire_matvec(A::SparseMatrixCSC{T,Int}, x::Vector{T}) where {T<:Posit}
    n = size(A, 2); y = Vector{T}(undef, n)
    rows = rowvals(A); vals = nonzeros(A)
    for j in 1:n
        q = Quire(T)
        for k in nzrange(A, j); fma_product!(q, vals[k], x[rows[k]]); end
        y[j] = T(q)
    end
    return y
end

function quire_matvec(A::SparseMatrixCSC{T,Int}, x::Vector{T}) where {T<:AnyBPosit}
    n = size(A, 2); y = Vector{T}(undef, n)
    rows = rowvals(A); vals = nonzeros(A)
    for j in 1:n
        q = BPosits.Quire()
        for k in nzrange(A, j); BPosits.add_product!(q, vals[k], x[rows[k]]); end
        y[j] = T(q)
    end
    return y
end

quire_matvec(A::SparseMatrixCSC{Float64,Int}, x::Vector{Float64}) = A * x



function load_matrix(mdname::String)
    # Do NOT wrap in Matrix(): densifying s3dkq4m2 (n=90449) is 65 GB.
    A0 = sparse(matrixdepot(mdname))

    eltype(A0) === Bool && error(
        "$mdname is a PATTERN-ONLY matrix (eltype Bool) with no numerical " *
        "values; use the valued copy in the other SuiteSparse group")
    A = eltype(A0) === Float64 ? A0 : SparseMatrixCSC{Float64,Int}(A0)
    issymmetric(A) || error("$mdname is not symmetric; CG requires SPD")
    if SYMSCALE
        d = 1 ./ sqrt.(diag(A))
        A = sparse(Diagonal(d) * A * Diagonal(d))
        A = (A + A') / 2
    end
    return A
end

function matrix_stats(A::SparseMatrixCSC{Float64})
    n  = size(A, 1)
    dn = filter(>(0), abs.(diag(A)))
    v  = filter(>(0), abs.(nonzeros(A)))
    spec = n <= SPECTRUM_MAX_N ? extrema(eigvals(Symmetric(Matrix(A)))) : (NaN, NaN)
    return (n=n, nnz=nnz(A),
            diag_min=minimum(dn), diag_max=maximum(dn), diag_ratio=maximum(dn)/minimum(dn),
            val_min=minimum(v), val_max=maximum(v), val_ratio=maximum(v)/minimum(v),
            lam_min=spec[1], lam_max=spec[2])
end


function run_ladder(name::String, A::SparseMatrixCSC{Float64}, x_true::Vector{Float64};
                    maxiter::Int=MAXITER, write_csv::Bool=true, tag_suffix::String="")
    n = size(A, 1)
    b = A * x_true
    diagA = diag(A)
    any(iszero, diagA) && error("$name has a zero diagonal; Jacobi PC undefined")

    x = zeros(Float64, n)
    r = copy(b)
    z = r ./ diagA
    p = copy(z)
    rz_d = naive_dot(r, z)

    relerrs = Dict{String,Vector{Float64}}()
    xc_same, xc_worst = Dict{String,Int}(), Dict{String,Float64}()
    for (_, _, tg) in LADDER, m in ("q", "n")
        relerrs["$(tg)$(m)"] = Float64[]; relerrs["$(tg)$(m)_bp"] = Float64[]
        xc_same["$(tg)$(m)"] = 0; xc_worst["$(tg)$(m)"] = 0.0
    end
    div_iter = Dict{String,Int}()   # first iter where naive and quire differ
  
    conv_iter = -1
    bnorm = sqrt(naive_dot(b, b))

    csv = write_csv ? open("results/ladder_$(name)$(tag_suffix).csv", "w") : nothing
    if write_csv
        hdr = ["iter", "pAp_d"]
        for (_, _, tg) in LADDER
            append!(hdr, ["$(tg)q","$(tg)n","$(tg)q_relerr","$(tg)n_relerr",
                          "$(tg)q_bp","$(tg)n_bp","$(tg)q_bp_relerr","$(tg)n_bp_relerr",
                          "$(tg)q_same","$(tg)n_same"])
        end
        push!(hdr, "residual")
        println(csv, join(hdr, ","))
    end

    for iter in 0:(maxiter - 1)
        Ap = A * p
        pAp_d = naive_dot(p, Ap)
        row = Any[iter, pAp_d]

        for (T, B, tg) in LADDER
            pT, ApT = T.(p), T.(Ap)
            pB, ApB = B.(p), B.(Ap)
            q_val = Float64(quire_dot(pT, ApT))
            n_val = Float64(naive_dot(pT, ApT))
            q_bp  = Float64(quire_dot(pB, ApB))
            n_bp  = Float64(naive_dot(pB, ApB))

            qe, ne   = relerr(q_val, pAp_d), relerr(n_val, pAp_d)
            qeb, neb = relerr(q_bp,  pAp_d), relerr(n_bp,  pAp_d)
            push!(relerrs["$(tg)q"], qe);     push!(relerrs["$(tg)n"], ne)
            push!(relerrs["$(tg)q_bp"], qeb); push!(relerrs["$(tg)n_bp"], neb)

            # Divergence onset: first iteration where quire and naive
            # relative errors separate by more than a factor of 2.
            if !haskey(div_iter, tg) && isfinite(qe) && isfinite(ne) && qe > 0 && ne / qe > 2
                div_iter[tg] = iter
            end

            qs, ns = same(q_val, q_bp), same(n_val, n_bp)
            xc_same["$(tg)q"] += qs; xc_same["$(tg)n"] += ns
            qs || (xc_worst["$(tg)q"] = max(xc_worst["$(tg)q"], relerr(q_bp, q_val)))
            ns || (xc_worst["$(tg)n"] = max(xc_worst["$(tg)n"], relerr(n_bp, n_val)))
            append!(row, [q_val, n_val, qe, ne, q_bp, n_bp, qeb, neb, qs, ns])
        end

        alpha = rz_d / pAp_d
        x = x .+ alpha .* p
        r = r .- alpha .* Ap
        z = r ./ diagA
        rz2 = naive_dot(r, z); beta = rz2 / rz_d; rz_d = rz2
        p = z .+ beta .* p
        rn = sqrt(naive_dot(r, r))
        conv_iter < 0 && rn / bnorm < RTOL && (conv_iter = iter)
        push!(row, rn)
        write_csv && println(csv, join(row, ","))
    end
    write_csv && close(csv)
    return (relerrs=relerrs, same=xc_same, worst=xc_worst, div_iter=div_iter,
            iters=maxiter, conv_iter=conv_iter,
            # window over which the ladder statistics are meaningful
            nstat = conv_iter < 0 ? maxiter : conv_iter + 1)
end


function pcg(A::AbstractMatrix{T}, b::Vector{T}, x_true::Vector{T};
             usepc::Bool, ipfun::Function, Aref::SparseMatrixCSC{Float64},
             mvfun::Function=plain_matvec,
             maxiter::Int=MAXITER, rtol::Float64=RTOL, keep_hist::Bool=true) where {T}
    n = size(A, 1)
    diagA_T = T[A[i, i] for i in 1:n]
    x = zeros(T, n)
    r = copy(b)
    z = usepc ? r ./ diagA_T : copy(r)
    p = copy(z)
    rz = ipfun(r, z)

    bnorm = sqrt(Float64(naive_dot(b, b)))
    xtd = Float64.(x_true)
    xtnorm = sqrt(naive_dot(xtd, xtd))
    hist = NamedTuple[]
    iters_to_converge = maxiter
    breakdown = "none"
    pd_prev, pApref_prev = Float64[], NaN

    errA_of(xv) = (e = Float64.(xv) .- xtd; q = naive_dot(e, Aref * e);
                   q < 0 ? NaN : sqrt(q))
    # ||x* - x_0||_A with x_0 = 0, i.e. ||x*||_A. Carson and Gergelits plot
    # ||x* - x_k||_A / ||x* - x_0||_A, which starts at 1 and is comparable
    # across matrices; errA alone carries the scale of x* and is not.
    xtA = (q0 = naive_dot(xtd, Aref * xtd); q0 <= 0 ? NaN : sqrt(q0))
    errA_rel_of(xv) = errA_of(xv) / xtA

    for k in 1:maxiter
        Ap  = mvfun(A, p)
        pAp = ipfun(p, Ap)

        pd      = Float64.(p)
        pAp_v   = Float64(pAp)
        pAp_dot = naive_dot(pd, Float64.(Ap))
        pAp_ref = naive_dot(pd, Aref * pd)
        ppd     = naive_dot(pd, pd)
        rayleigh     = iszero(ppd) ? NaN : pAp_v / ppd
        rayleigh_ref = iszero(ppd) ? NaN : pAp_ref / ppd
        conjloss = if isempty(pd_prev) || !(pAp_ref > 0 && pApref_prev > 0)
            NaN
        else
            abs(naive_dot(pd, Aref * pd_prev)) / sqrt(pAp_ref * pApref_prev)
        end

        if !isfinite(pAp_v) || pAp_v <= 0
            breakdown = !isfinite(pAp_v) ? "pAp_nonfinite" : "pAp_nonpositive"
            keep_hist && push!(hist, (iter=k, residual=NaN, error=NaN, pAp=pAp_v,
                pAp_dot=pAp_dot, pAp_ref=pAp_ref,
                relerr_dot=relerr(pAp_v, pAp_dot), relerr_stp=relerr(pAp_v, pAp_ref),
                rayleigh=rayleigh, rayleigh_ref=rayleigh_ref,
                conjloss=conjloss, posdef=false, errA=errA_of(x),
                errA_rel=errA_rel_of(x)))
            iters_to_converge = k
            break
        end

        alpha = rz / pAp
        x = x .+ alpha .* p
        r = r .- alpha .* Ap

        resnorm = sqrt(Float64(naive_dot(r, r))) / bnorm
        errnorm = sqrt(sum((Float64.(x) .- xtd) .^ 2))
        keep_hist && push!(hist, (iter=k, residual=resnorm, error=errnorm, pAp=pAp_v,
            pAp_dot=pAp_dot, pAp_ref=pAp_ref,
            relerr_dot=relerr(pAp_v, pAp_dot), relerr_stp=relerr(pAp_v, pAp_ref),
            rayleigh=rayleigh, rayleigh_ref=rayleigh_ref,
            conjloss=conjloss, posdef=true, errA=errA_of(x),
            errA_rel=errA_rel_of(x)))

        if resnorm < rtol
            iters_to_converge = k
            break
        end
        pd_prev, pApref_prev = pd, pAp_ref
        z  = usepc ? r ./ diagA_T : copy(r)
        rz2 = ipfun(r, z); beta = rz2 / rz; rz = rz2
        p = z .+ beta .* p
    end

    xd = Float64.(x)
    return (x=x, iters=iters_to_converge, hist=hist, breakdown=breakdown,
            rel_error = sqrt(naive_dot(xd .- xtd, xd .- xtd)) / xtnorm,
            errA = errA_of(x), errA_rel = errA_rel_of(x))
end

const VARIANTS = [
    ("Float64",   Float64,             naive_dot, "naive"),
    ("Posit8_2",  Posit{8,2,UInt8},    naive_dot, "naive"),
    ("Posit8_2",  Posit{8,2,UInt8},    quire_dot, "quire"),
    ("BPosit8",   BPosit8,             naive_dot, "naive"),
    ("BPosit8",   BPosit8,             quire_dot, "quire"),
    ("Posit16_2", Posit{16,2,UInt16},  naive_dot, "naive"),
    ("Posit16_2", Posit{16,2,UInt16},  quire_dot, "quire"),
    ("BPosit16",  BPosit16,            naive_dot, "naive"),
    ("BPosit16",  BPosit16,            quire_dot, "quire"),
    ("Posit32_2", Posit{32,2,UInt32},  naive_dot, "naive"),
    ("Posit32_2", Posit{32,2,UInt32},  quire_dot, "quire"),
    ("BPosit32",  BPosit32,            naive_dot, "naive"),
    ("BPosit32",  BPosit32,            quire_dot, "quire"),
]



const SEED_VARIANTS = [
    ("Float64",   Float64,             naive_dot, "naive"),
    ("Posit32_2", Posit{32,2,UInt32},  naive_dot, "naive"),
    ("Posit32_2", Posit{32,2,UInt32},  quire_dot, "quire"),
]

const HIST_COLS = (:iter, :residual, :error, :pAp, :pAp_dot, :pAp_ref,
                   :relerr_dot, :relerr_stp, :rayleigh, :rayleigh_ref,
                   :conjloss, :posdef, :errA, :errA_rel)

function run_driven(name::String, A::SparseMatrixCSC{Float64}, x_true::Vector{Float64},
                    st, io::IO)
    n = size(A, 1)
    mkpath("results/$name")
    for (label, T, ipfun, ipname) in VARIANTS
        A_T = T == Float64 ? A : SparseMatrixCSC{T,Int}(A)
        xt = T.(x_true)
        bT = A_T * xt
        for usepc in (true, false)
            res = pcg(A_T, bT, xt; usepc=usepc, ipfun=ipfun, Aref=A)
            pcname = usepc ? "jacobi" : "none"
            open("results/$name/pcg_$(label)_$(ipname)_$(pcname).csv", "w") do f
                println(f, join(HIST_COLS, ","))
                for h in res.hist
                    println(f, join((getfield(h, c) for c in HIST_COLS), ","))
                end
            end
            fin(f) = isempty(res.hist) ? NaN : getfield(res.hist[end], f)
            finite_of(f) = [getfield(h, f) for h in res.hist if isfinite(getfield(h, f))]
            maxof(f) = (v = finite_of(f); isempty(v) ? NaN : maximum(v))
            eA = [h.errA for h in res.hist]
            rises = [k for k in 2:length(eA)
                     if isfinite(eA[k]) && isfinite(eA[k-1]) && eA[k] > eA[k-1]]
            bad = findfirst(h -> !h.posdef, res.hist)
            rq = finite_of(:rayleigh)
            in_spec = isnan(st.lam_min) ? "NA" : isempty(rq) ? "false" :
                      string(minimum(rq) >= st.lam_min*(1-1e-8) &&
                             maximum(rq) <= st.lam_max*(1+1e-8))
            println(io, "$name,$label,$ipname,$pcname,$(res.iters),$(res.breakdown)," *
                        "$(fin(:residual)),$(res.rel_error),$(res.errA),$(res.errA_rel)," *
                        "$(length(rises)),$(isempty(rises) ? "" : first(rises))," *
                        "$(maxof(:relerr_dot)),$(maxof(:conjloss))," *
                        "$(bad === nothing ? "" : bad),$in_spec")
            flush(io)
        end
    end
end

# =================================================================
# Stage 3: seed study over random GROUND TRUTHS (see GROUND TRUTH note)
# =================================================================
function run_seed_study(name::String, A::SparseMatrixCSC{Float64}, st,
                        seed_io::IO, agg_io::IO)
    n = size(A, 1)
    do_driven = st.nnz <= SEED_DRIVEN_MAX_NNZ
    div_by_tag = Dict(tg => Int[] for (_, _, tg) in LADDER)
    ndet = Dict(tg => 0 for (_, _, tg) in LADDER)
    drv = Dict{String,Vector{Float64}}()
    brk = Dict{String,Int}()

    for seed in 1:NSEEDS
        rng = MersenneTwister(seed)
        x_true = randn(rng, n)          # random ground truth, b = A*x
        lad = run_ladder(name, A, x_true; write_csv=false)

        for (_, _, tg) in LADDER
            di = get(lad.div_iter, tg, -1)
            di >= 0 && (push!(div_by_tag[tg], di); ndet[tg] += 1)
            println(seed_io, "$name,$seed,$tg,ladder,$(di < 0 ? "" : di),,,")
        end

        if do_driven
            for (label, T, ipfun, ipname) in SEED_VARIANTS
                A_T = T == Float64 ? A : SparseMatrixCSC{T,Int}(A)
                xt = T.(x_true); bT = A_T * xt
                res = pcg(A_T, bT, xt; usepc=true, ipfun=ipfun, Aref=A, keep_hist=false)
                k = "$(label)_$(ipname)"
                push!(get!(drv, k, Float64[]), res.rel_error)
                res.breakdown != "none" && (brk[k] = get(brk, k, 0) + 1)
                println(seed_io, "$name,$seed,$label,$ipname,$(res.iters)," *
                                 "$(res.breakdown),$(res.rel_error),$(res.errA),$(res.errA_rel)")
            end
        end
        flush(seed_io)
    end

    for (_, _, tg) in LADDER
        d = div_by_tag[tg]
        println(agg_io, "$name,$tg,ladder_divergence,$NSEEDS,$(ndet[tg])," *
                        "$(100*ndet[tg]/NSEEDS),$(isempty(d) ? "" : mean(d))," *
                        "$(length(d) > 1 ? std(d) : "")," *
                        "$(isempty(d) ? "" : minimum(d)),$(isempty(d) ? "" : maximum(d))")
    end
    for (k, v) in sort(collect(drv), by=first)
        println(agg_io, "$name,$k,driven_rel_error,$NSEEDS,$(get(brk,k,0))," *
                        "$(100*get(brk,k,0)/NSEEDS),$(mean(v)),$(std(v))," *
                        "$(minimum(v)),$(maximum(v))")
    end
    flush(agg_io)
    return do_driven
end


# =================================================================
# Stage 4: is the MATVEC the throttle?

const MV_GRID = [("plain", plain_matvec), ("quire", quire_matvec)]
# Posit16_2 dropped: it produces relative solution errors of 1e0 to 1e12 on
# every matrix here, so no 16-bit cell of the factorial is interpretable.
const S4_TYPES = [("Posit32_2", Posit{32,2,UInt32}), ("BPosit32", BPosit32)]

function run_stage4(name::String, A::SparseMatrixCSC{Float64}, io::IO, seeded::Bool)
    n = size(A, 1)
    ns = nnz(A) > S4_LARGE_NNZ ? NSEEDS_S4_LARGE : NSEEDS
    seeds = seeded ? (1:ns) : (1:1)
    for (label, T) in S4_TYPES
        A_T = SparseMatrixCSC{T,Int}(A)
        for seed in seeds
            x_true = seeded ? randn(MersenneTwister(seed), n) : ones(n)
            xt = T.(x_true)
            bT = T.(A * x_true)          # common RHS, see note above
            for (dotname, ipfun) in (("naive", naive_dot), ("quire", quire_dot)),
                (mvname, mvfun) in MV_GRID
                res = pcg(A_T, bT, xt; usepc=true, ipfun=ipfun, Aref=A,
                          mvfun=mvfun, keep_hist=false)
                println(io, "$name,$label,$seed,$dotname,$mvname,$(res.iters)," *
                            "$(res.breakdown),$(res.rel_error),$(res.errA),$(res.errA_rel)")
            end
        end
        flush(io)
    end
end

# =================================================================
# Driver
# =================================================================
mkpath("results")
stats_io = open("results/matrix_stats.csv", "w")
println(stats_io, "matrix,mdname,n,nnz,diag_min,diag_max,diag_ratio,val_min,val_max," *
                  "val_ratio,lam_min,lam_max,conv_iter,stat_window,driven,seed_study," *
                  "seed_driven,ladder_secs,driven_secs,seed_secs")
xc_io = open("results/crosscheck.csv", "w")
println(xc_io, "matrix,width,mode,iters,same_value,differ,worst_reldiff")
lad_io = open("results/ladder_summary.csv", "w")
println(lad_io, "matrix,width,impl,max_relerr_naive,max_relerr_quire,max_ratio," *
                "median_relerr_naive,median_relerr_quire,median_ratio,divergence_iter," *
                "conv_iter,stat_window")
sum_io = open("results/pcg_summary.csv", "w")
println(sum_io, "matrix,type,ipmode,preconditioner,iters,breakdown,final_residual," *
                "final_rel_error,final_errA,final_errA_rel," *
                "n_errA_increases,first_errA_increase," *
                "max_relerr_dot,max_conjloss,first_nonposdef_iter,rayleigh_in_spectrum")
seed_io = open("results/seed_runs.csv", "w")
println(seed_io, "matrix,seed,label,mode,iters,breakdown,rel_error,errA,errA_rel")
agg_io = open("results/seed_summary.csv", "w")
println(agg_io, "matrix,key,metric,n_seeds,n_detected,pct,mean,std,min,max")
s4_io = open("results/quired_matvec.csv", "w")
println(s4_io, "matrix,type,seed,dot,matvec,iters,breakdown,rel_error,errA,errA_rel")

@printf("%-10s %7s %9s %11s %9s %7s %9s %9s %9s\n", "matrix", "n", "nnz",
        "diag_ratio", "conv_iter", "window", "ladder_s", "driven_s", "seeds_s")
flush(stdout)

for (name, mdname) in MATRICES
    A = load_matrix(mdname)
    st = matrix_stats(A)
    x_ones = ones(Float64, st.n)

    t_lad = @elapsed lad = run_ladder(name, A, x_ones)
    for (_, _, tg) in LADDER, (m, mn) in (("q","quire"), ("n","naive"))
        k = "$(tg)$(m)"; s = lad.same[k]
        println(xc_io, "$name,$tg,$mn,$(lad.iters),$s,$(lad.iters-s),$(lad.worst[k])")
    end
    for (_, _, tg) in LADDER, (impl, sfx) in (("UniversalNumbers",""), ("BPosits","_bp"))
        w  = lad.nstat
        nv = filter(isfinite, lad.relerrs["$(tg)n$(sfx)"][1:w])
        qv = filter(isfinite, lad.relerrs["$(tg)q$(sfx)"][1:w])
        mn = isempty(nv) ? NaN : maximum(nv); mq = isempty(qv) ? NaN : maximum(qv)
        dn = isempty(nv) ? NaN : median(nv);  dq = isempty(qv) ? NaN : median(qv)
        println(lad_io, "$name,$tg,$impl,$mn,$mq,$(mn/mq),$dn,$dq,$(dn/dq)," *
                        "$(get(lad.div_iter, tg, "")),$(lad.conv_iter),$w")
    end
    flush(xc_io); flush(lad_io)

    driven = st.nnz <= DRIVEN_MAX_NNZ
    t_drv = driven ? @elapsed(run_driven(name, A, x_ones, st, sum_io)) : NaN

    seeds = st.n <= SEED_MAX_N
    seed_driven = false
    t_seed = NaN
    if seeds
        t_seed = @elapsed (seed_driven = run_seed_study(name, A, st, seed_io, agg_io))
    end

    # Stage 4: seeded where the seed budget allows, single-shot otherwise
    st.nnz <= DRIVEN_MAX_NNZ && run_stage4(name, A, s4_io, st.nnz <= SEED_DRIVEN_MAX_NNZ)

    println(stats_io, "$name,$mdname,$(st.n),$(st.nnz),$(st.diag_min),$(st.diag_max)," *
                      "$(st.diag_ratio),$(st.val_min),$(st.val_max),$(st.val_ratio)," *
                      "$(st.lam_min),$(st.lam_max),$(lad.conv_iter),$(lad.nstat)," *
                      "$driven,$seeds,$seed_driven," *
                      "$t_lad,$t_drv,$t_seed")
    flush(stats_io)
    @printf("%-10s %7d %9d %11.3e %9s %7d %9.1f %9s %9s\n", name, st.n, st.nnz,
            st.diag_ratio, lad.conv_iter < 0 ? ">300" : string(lad.conv_iter),
            lad.nstat, t_lad, driven ? @sprintf("%.1f", t_drv) : "skip",
            seeds ? @sprintf("%.1f", t_seed) : "skip")
    flush(stdout)
    A = nothing; GC.gc()
end

close(stats_io); close(xc_io); close(lad_io); close(sum_io); close(seed_io)
close(agg_io); close(s4_io)
println("\nDone. results/: matrix_stats.csv ladder_summary.csv crosscheck.csv " *
        "pcg_summary.csv seed_runs.csv seed_summary.csv ladder_<m>.csv <m>/pcg_*.csv")
