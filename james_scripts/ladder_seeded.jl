using UniversalNumbers
using SparseArrays
using LinearAlgebra
using MatrixDepot
using Random
using Statistics
using Printf


const MAXITER = 300
const RTOL    = 1e-10
const NSEEDS  = 50
const MAX_N   = 90_449

const WIDTHS = [("p8",  Posit{8,2,UInt8}),
                ("p16", Posit{16,2,UInt16}),
                ("p32", Posit{32,2,UInt32})]

const MATRICES = [
    ("bcsstk03", "HB/bcsstk03"),        ("nos2",     "HB/nos2"),
    ("mhd4800b", "Bai/mhd4800b"),       ("bcsstk14", "HB/bcsstk14"),
    ("sts4098",  "Cannizzo/sts4098"),   ("nasa4704", "Nasa/nasa4704"),
    ("bodyy4",   "Pothen/bodyy4"),      ("bcsstk38", "Boeing/bcsstk38"),
    ("bcsstk37", "Boeing/bcsstk37"),    ("bcsstk36", "Boeing/bcsstk36"),
    ("nasasrb",  "Nasa/nasasrb"),       ("s3dkt3m2", "GHS_psdef/s3dkt3m2"),
    ("s3dkq4m2", "GHS_psdef/s3dkq4m2"),
]

function naive_dot(a::AbstractVector, b::AbstractVector)
    s = zero(eltype(a))
    for i in eachindex(a, b)
        s += a[i] * b[i]
    end
    return s
end
quire_dot(a::AbstractVector{T}, b::AbstractVector{T}) where {T<:Posit} = fdp(a, b)

relerr(v, ref) = iszero(ref) ? (iszero(v) ? 0.0 : Inf) : abs(v - ref) / abs(ref)
safe(f, v) = (w = filter(isfinite, v); isempty(w) ? NaN : f(w))

"""
One double-driven ladder pass. Returns, per width and per reduction
(pAp and rz), the vectors of relative errors over the pre-convergence
window, plus the convergence iteration.
"""
function ladder_pass(A::SparseMatrixCSC{Float64}, x_true::Vector{Float64})
    n = size(A, 1)
    b = A * x_true
    d = diag(A)
    x = zeros(n); r = copy(b); z = r ./ d; p = copy(z)
    rz_d = naive_dot(r, z)
    bn = sqrt(naive_dot(b, b))

    E = Dict{String,Vector{Float64}}()
    for (tag, _) in WIDTHS, red in ("pAp", "rz"), mode in ("n", "q")
        E["$(tag)_$(red)_$(mode)"] = Float64[]
    end
    E["errA_rel"] = Float64[]
    q0 = naive_dot(x_true, A * x_true)
    nA0 = q0 <= 0 ? NaN : sqrt(q0)
    conv = -1

    for iter in 0:(MAXITER - 1)
        Ap = A * p
        pAp_d = naive_dot(p, Ap)
        ed = x .- x_true
        qe = naive_dot(ed, A * ed)
        push!(E["errA_rel"], qe < 0 ? NaN : sqrt(qe) / nA0)

        for (tag, T) in WIDTHS
            pT, ApT = T.(p), T.(Ap)
            push!(E["$(tag)_pAp_n"], relerr(Float64(naive_dot(pT, ApT)), pAp_d))
            push!(E["$(tag)_pAp_q"], relerr(Float64(quire_dot(pT, ApT)), pAp_d))
            rT, zT = T.(r), T.(z)
            push!(E["$(tag)_rz_n"], relerr(Float64(naive_dot(rT, zT)), rz_d))
            push!(E["$(tag)_rz_q"], relerr(Float64(quire_dot(rT, zT)), rz_d))
        end

        al = rz_d / pAp_d
        x = x .+ al .* p
        r = r .- al .* Ap
        z = r ./ d
        rz2 = naive_dot(r, z)
        p = z .+ (rz2 / rz_d) .* p
        rz_d = rz2
        conv < 0 && sqrt(naive_dot(r, r)) / bn < RTOL && (conv = iter)
    end

    w = conv < 0 ? MAXITER : conv + 1       # statistics window
    for k in keys(E); E[k] = E[k][1:w]; end
    return E, conv, w
end

mkpath("results")
runs = open("results/ladder_seeded_runs.csv", "w")
println(runs, "matrix,seed,width,reduction,conv_iter,window," *
              "max_naive,max_quire,max_ratio,med_naive,med_quire,med_ratio," *
              "errA_rel_final")
agg = open("results/ladder_seeded_summary.csv", "w")
println(agg, "matrix,width,reduction,n_seeds,mean_max_ratio,std_max_ratio," *
             "min_max_ratio,max_max_ratio,median_max_ratio," *
             "mean_med_ratio,std_med_ratio,mean_conv_iter")

@printf("%-10s %-5s %-5s %10s %10s %10s %10s\n",
        "matrix", "width", "red", "mean_ratio", "std", "min", "max")
flush(stdout)

for (name, mdname) in MATRICES
    A0 = sparse(matrixdepot(mdname))
    eltype(A0) === Bool && continue
    A = eltype(A0) === Float64 ? A0 : SparseMatrixCSC{Float64,Int}(A0)
    n = size(A, 1)
    if n > MAX_N
        @printf("%-10s n=%d exceeds MAX_N=%d, skipped\n", name, n, MAX_N)
        flush(stdout); A = nothing; A0 = nothing; GC.gc(); continue
    end

    acc = Dict{Tuple{String,String},Vector{Float64}}()   # (width,red) -> ratios
    accm = Dict{Tuple{String,String},Vector{Float64}}()
    convs = Int[]

    for seed in 1:NSEEDS
        x_true = randn(MersenneTwister(seed), n)
        E, conv, w = ladder_pass(A, x_true)
        push!(convs, conv < 0 ? MAXITER : conv)
        eAf = isempty(E["errA_rel"]) ? NaN : last(E["errA_rel"])
        for (tag, _) in WIDTHS, red in ("pAp", "rz")
            mn = safe(maximum, E["$(tag)_$(red)_n"])
            mq = safe(maximum, E["$(tag)_$(red)_q"])
            dn = safe(median,  E["$(tag)_$(red)_n"])
            dq = safe(median,  E["$(tag)_$(red)_q"])
            rmax, rmed = mn / mq, dn / dq
            isfinite(rmax) && push!(get!(acc, (tag, red), Float64[]), rmax)
            isfinite(rmed) && push!(get!(accm, (tag, red), Float64[]), rmed)
            println(runs, "$name,$seed,$tag,$red,$conv,$w,$mn,$mq,$rmax,$dn,$dq,$rmed,$eAf")
        end
    end
    flush(runs)

    for (tag, _) in WIDTHS, red in ("pAp", "rz")
        v = get(acc, (tag, red), Float64[]); isempty(v) && continue
        vm = get(accm, (tag, red), Float64[])
        println(agg, "$name,$tag,$red,$(length(v)),$(mean(v)),$(std(v))," *
                     "$(minimum(v)),$(maximum(v)),$(median(v))," *
                     "$(isempty(vm) ? NaN : mean(vm)),$(length(vm)>1 ? std(vm) : NaN)," *
                     "$(mean(convs))")
        @printf("%-10s %-5s %-5s %10.4g %10.4g %10.4g %10.4g\n",
                name, tag, red, mean(v), std(v), minimum(v), maximum(v))
        flush(stdout)
    end
    flush(agg)
    A = nothing; A0 = nothing; GC.gc()
end

close(runs); close(agg)
println("\nDone -> results/ladder_seeded_runs.csv, results/ladder_seeded_summary.csv")
