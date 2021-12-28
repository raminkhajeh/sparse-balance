module v
using JLD
using FFTW
using Trapz
using Printf
using PyPlot
using Random
using Statistics
using Parameters
using SparseArrays
using ProgressMeter
using Distributions
using LinearAlgebra
using DelimitedFiles

@with_kw mutable struct Params
   prefix::String = "D"
   record::Bool = true        # record trajectories if true
   compACF::Bool = true       # compute autocorrelation function if true

   N::Int = 1000              # number of neurons
   K::Int = 1000              # in-degree
   T::Int = 2000              # simulation time
   warmUp::Int = 10           # warm up time
   ncorr::Int = 40            # autocorrelation time
   dt::Float64 = 0.02         # time step
   τ::Float64 = 1.0           # time constant 
   rcNsamp::Int = 1000;       @assert rcNsamp <= N;   # recorded sample
   rcDt::Float64 = 0.02;      @assert rcDt >= dt;     # recorded time step
   Φname::String = "rectTanh" # response function name
   λ::Float64 = 1.0           # powerlaw exponent
   Jdist::String = "gamma"    # connectivity distribution name
   ρ_e::Float64 = 0.          # fraction of excitatory neurons

   J_ee::Float64 = 0.         # Je→e (sets the mean coupling)
   J_ei::Float64 = 0.         # Ji→e
   J_ie::Float64 = 0.         # Je→i
   J_ii::Float64 = 2.0        # Ji→i
   scaleJ::Float64 = -0.5     # scaling for mean coupling 

   g_ee::Float64 = 0.         # ge→e (sets the std of coupling)
   g_ei::Float64 = 0.         # gi→e 
   g_ie::Float64 = 0.         # ge→i
   g_ii::Float64 = 2.0        # gi→i
   scaleg::Float64 = -0.25    # scaling for std coupling

   I_e::Float64 = 0.          # input current to e
   I_i::Float64 = 1.0         # input current to i
   scaleI::Float64 = 0.0      # scaling for input current

   seedJ::Int = 0             # network connectivity seed
   seedIC::Int = 0            # initial condition seed
   progress::Bool = true      # show progress if true

end
dp = Params()

#################################################
function runsim(;p = dp)
   startTime = time()
   function rectTanh(x)
      if x > 0
         return tanh(x)
      else
         return 0.0
      end
   end
   function rectLin(x)
      if x > 0
         return x
      else
         return 0.0
      end
   end
   function Heaviside(x)
      if x > 0
         return 1.0
      else
         return 0.0
      end
   end
   function rectPower(x)
      if x > 0
         return x^p.λ
      else
         return 0.0
      end
   end
   diluteJ = (p.N != p.K)

   N_i = round(Int, p.N - p.N * p.ρ_e)
   N_e = p.N - N_i
   K_e = round(Int, p.K * p.ρ_e)
   K_i = round(Int, p.K - p.K * p.ρ_e)

   I_e = p.I_e * (K_e ^ p.scaleI)
   I_i = p.I_i * (K_i ^ p.scaleI)

   tstep = floor(Integer, p.T / p.dt)
   tstep_warmUp = floor(Integer, p.warmUp * p.τ / p.dt)
   tstep_tot = tstep + tstep_warmUp
   tcorr_step = floor(Integer, p.ncorr * p.τ / p.rcDt)
   tcorr = 0 : p.rcDt : p.ncorr * p.τ - p.rcDt
   τ⁻¹dt = p.dt / p.τ

   ΦDict = Dict("rectTanh" => rectTanh, "rectLin" => rectLin, 
                "rectPower" => rectPower, "Heaviside" => Heaviside)
   Φ = ΦDict[p.Φname]

   rcEvery = floor(Int, p.rcDt/p.dt)
   rctstep = floor(Int, tstep/rcEvery)

   if (p.ρ_e != 0)
      rcN_e        = Int(p.ρ_e * p.rcNsamp)
      rcNsampE     = 1:rcN_e
      rcNsampI     = (rcN_e + 1): p.rcNsamp
      rcNsampRange = vcat(rcNsampE,rcNsampI)
   else
      rcNsampRange = 1:p.rcNsamp
   end

   D = Dict()
   F = Dict()

   if (diluteJ)
      if !iszero(p.ρ_e)
         printstyled("partially-connected EI","\n";color=:green)
         J = randJ_EI_PC(p = p)
      else
         printstyled("partially-connected pure-I","\n";color=:green)
         J = randJ_I_PC(p = p)
      end
   else
      if !iszero(p.ρ_e)
         printstyled("fully-connected EI","\n";color=:green)
         J = randJ_EI_FC(p = p)
      else
         printstyled("fully-connected pure-I","\n";color=:green)
         J = randJ_I_FC(p = p)
      end
   end
   rcCount = 0

   x       = p.g_ii * randn(MersenneTwister(p.seedIC), p.N)
   ϕ       = zeros(Float64,p.N)
   I       = zeros(Float64,p.N)
   percTon = zeros(Int,p.N)

   m   = zeros(Float64,rctstep)
   x̄   = zeros(Float64,rctstep)
   η̄   = zeros(Float64,rctstep)
   n   = zeros(Int,rctstep)
   μ   = zeros(Float64,rctstep)

   if (p.record || p.compACF)
      xs = zeros(Float64,rctstep, p.rcNsamp)
      ϕs = zeros(Float64,rctstep, p.rcNsamp)
      ηs = zeros(Float64,rctstep, p.rcNsamp)
   end

   if !iszero(p.ρ_e)
      mE = zeros(Float64,rctstep)
      mI = zeros(Float64,rctstep)
      nE = zeros(Int,rctstep) 
      nI = zeros(Int,rctstep)
   end

   tmx  = zeros(Float64,p.N)
   tmη  = zeros(Float64,p.N)
   tmx2 = zeros(Float64,p.N)
   tmη2 = zeros(Float64,p.N)

   I[1:N_e] .= I_e
   I[(N_e + 1):p.N] .= I_i
   ϕ = Φ.(x)
   iOn = findall(ϕ .> 0)

   ##########################
   t = 0.0
   printParams(J, I_i, p = p)
   prog = Progress(tstep_tot, 0.5)
   for ts = 1:tstep_tot
      p.progress && ProgressMeter.update!(prog,ts);
      t₀ = ts - tstep_warmUp
      #######################
      η = J * ϕ
      x = x + τ⁻¹dt * (-x + η + I)
      ϕ = Φ.(x)
      #######################
      if (t₀ > 0) && (t₀ % rcEvery == 0)
         iOn = findall(x .> 0)
         percTon[iOn] .+= 1
         rcCount += 1

         n[rcCount]   = length(iOn)
         η̄[rcCount]   = mean(η)
         x̄[rcCount]   = mean(x)
         m[rcCount]   = mean(ϕ)
         μ[rcCount]   = mean(ϕ[iOn])
         tmx  .+= x
         tmη  .+= η
         tmx2 .+= x.^2
         tmη2 .+= η.^2

         if !iszero(p.ρ_e)
            mE[rcCount] = mean(ϕ[1:N_e])
            mI[rcCount] = mean(ϕ[N_e+1:end])
            nE[rcCount] = count(iOn .<= N_e)
            nI[rcCount] = count(iOn .> N_e)
         end

         if (p.record || p.compACF)
            xs[rcCount,:] = x[rcNsampRange]
            ϕs[rcCount,:] = ϕ[rcNsampRange]
            ηs[rcCount,:] = η[rcNsampRange]
         end
      end
   end

   #####################
   D[:m] = mean(m)                     # [⟨ϕ⟩]
   D[:x̄] = mean(x̄)                     # [⟨x⟩]
   D[:η̄] = mean(η̄)                     # [⟨η⟩]

   ####################
   D[:Mtm2x] = mean((tmx/rcCount).^2)  # [⟨x⟩²]
   D[:Mtm2η] = mean((tmη/rcCount).^2)  # [⟨η⟩²]

   D[:Mtmx2] = mean((tmx2/rcCount))    # [⟨x²⟩]
   D[:Mtmη2] = mean((tmη2/rcCount))    # [⟨η²⟩]

   D[:tvx] = D[:Mtmx2] - D[:Mtm2x]     # temporal var of x
   D[:tvη] = D[:Mtmη2] - D[:Mtm2η]     # temporal var of η

   D[:qvx] = D[:Mtm2x] - D[:x̄]^2       # quenched var of x
   D[:qvη] = D[:Mtm2η] - D[:η̄]^2       # quenched var of η

   #####################
   D[:f̄], D[:Δf²] = mean(n/p.N), var(n/p.N)
   D[:μ̄], D[:Δμ²] = mean(μ), var(μ)
   D[:percTon] = percTon/rctstep

   #####################
   D[:α] = mean(m) * sqrt(p.K) * (p.J_ii / p.g_ii)^2

   #####################
   if !iszero(p.ρ_e)
      D[:mE] = mE
      D[:mI] = mI
      D[:fE] = nE / N_e
      D[:fI] = nI / N_i
   end

   #####################
   if (p.compACF)
      Rδη = comp_ACF(ηs .- mean(ηs,dims=1), tcorr_step);D[:Rδη]=Rδη
      Rδϕ = comp_ACF(ϕs .- mean(ϕs,dims=1), tcorr_step);D[:Rδϕ]=Rδϕ
      Rδx = comp_ACF(xs .- mean(xs,dims=1), tcorr_step);D[:Rδx]=Rδx

      D[:β] = comp_β(tcorr, Rδη; p = p)
   end

   if (p.record) 
      F[:ηs] = ηs
      F[:xs] = xs
      F[:ϕs] = ϕs
   else
      F[:xs] = xs[:,1:10]
      F[:ηs] = ηs[:,1:10]
   end

   D[:tcorr] = tcorr
   D[:δtime] = time() - startTime
   D[:rcNsampRange] = rcNsampRange
   !iszero(p.ρ_e) && (D[:rcNsampE]=rcNsampE; D[:rcNsampI]=rcNsampI)

   return D, F
end

#################################################
function randJ_generic(EorI, J̄, σJ, n, m; p = dp)
   if p.Jdist == "normal"
      J = EorI * J̄ .+ σJ * randn(n, m)
   elseif p.Jdist == "gamma"
      J = EorI * genGamma(J̄, σJ, n, m)
   elseif p.Jdist == "lognormal"
      J = EorI * genLognormal(J̄, σJ, n, m)
   elseif p.Jdist == "bernoulli"
      J = EorI * genBernoulli(J̄, n, m)
   end
   return J
end

#################################################
function randJ_I_FC(;p = dp)
   Random.seed!(p.seedJ)
   J̄ = p.J_ii * (p.K ^ p.scaleJ)
   σJ = p.g_ii * (p.K ^ p.scaleg)
   J = randJ_generic(-1, J̄, σJ, p.N, p.N, p = p)
   return J
end

#################################################
function randJ_I_PC(;p = dp)
   Random.seed!(p.seedJ)
   spMask = sprand(Bool, p.N, p.N, p.K/p.N)
   J_FC = randJ_I_FC(p = p)
   J = spMask .* J_FC
   return J
end

#################################################
function randJ_EI_FC(;p = dp)
   Random.seed!(p.seedJ)
   Ki = round(Int, p.K - p.K * p.ρ_e)
   Ke = p.K - Ki
   Ni = round(Int, p.N - p.N * p.ρ_e)
   Ne = p.N - Ni

   J̄ee = p.J_ee * (Ke ^ p.scaleJ)
   σJee = p.g_ee * (Ke ^ p.scaleg)
   Jee = randJ_generic(1,J̄ee,σJee,Ne,Ne,p=p)

   J̄ei = p.J_ei * (Ki ^ p.scaleJ)
   σJei = p.g_ei * (Ki ^ p.scaleg)
   Jei = randJ_generic(-1,J̄ei,σJei,Ne,Ni,p=p)

   J̄ie = p.J_ie * (Ke ^ p.scaleJ)
   σJie = p.g_ie * (Ke ^ p.scaleg)
   Jie = randJ_generic(1,J̄ie,σJie,Ni,Ne,p=p)

   J̄ii = p.J_ii * (Ki ^ p.scaleJ)
   σJii = p.g_ii * (Ki ^ p.scaleg)
   Jii = randJ_generic(-1,J̄ii,σJii,Ni,Ni,p=p)

   JfromE = vcat(Jee,Jie)
   JfromI = vcat(Jei,Jii)
   J = hcat(JfromE, JfromI)
   return J
end

#################################################
function randJ_EI_PC(;p = dp)
   Random.seed!(p.seedJ)
   spMask = sprand(Bool, p.N, p.N, p.K/p.N)
   J_FC = randJ_EI_FC(p = p)
   J = spMask .* J_FC
   return J
end

#################################################
function genBernoulli(μ, dim...)
   d = Bernoulli(μ)
   return rand(d, dim...)
end

#################################################
function genLognormal(μ, σ, dim...)
   M = log(μ^2 / sqrt(μ^2 + σ^2))
   Σ = sqrt(log((σ^2 / μ^2) + 1))
   d = LogNormal(M, Σ)
   return rand(d, dim...)
end

#################################################
function genGamma(μ, σ, dim...)
   θ = σ^2 / μ
   k = (μ / σ)^2
   d = Gamma(k, θ)
   return rand(d, dim...)
end

#################################################
uniform(a,b,n) = rand(n)*(b - a) .+ a
unroll(x) = reshape(x, length(x))
nanmean(x) = mean(filter(!isnan,x))
nanmean(x,dim) = dropdims(mapslices(nanmean,x,dims=dim),dims=dim)
nanstd(x) = std(filter(!isnan,x))
nanstd(x,dim) = dropdims(mapslices(nanstd,x,dims=dim),dims=dim)

#################################################
function myautocorr(x, maxLag; demean=true, normed=true)
   L = length(x)
   x = demean ? x .-= mean(x) : x
   yApn = zeros(2*L)
   yApn[1:L] .= x
   corr = irfft(abs.(rfft(yApn)).^2, length(yApn))
   B = floor(Int, length(corr)/2)
   corr = corr[1:B]
   lengths = B:-1:1
   corr = (corr ./ lengths)[1:maxLag]
   if var(x) == 0.0
      printstyled("Traj is flat.\n",color=:red)
      return ones(length(corr))
   else
      return (normed) ? corr / corr[1] : corr
   end
end

#################################################
function comp_ACF(y, maxLag)
   N = size(y, 2)
   nonzeroN = findall(!iszero, dropdims(var(y,dims=1),dims=1))
   ac = zeros(maxLag)
   for i in nonzeroN
      ac .+= myautocorr(y[:,i], maxLag, demean=false, normed=false)
   end
   ac /= length(nonzeroN)
   return ac
end

#################################################
function comp_β(τ, Rδη; p = dp)
   β = Rδη[1] ./ trapz(τ, Rδη)
   return β
end

#################################################
function printParams(J, I_i;p = dp)
   runparams = string(
            "N = ",       @sprintf("%d",p.N),
         " | K = ",       @sprintf("%d",p.K),
         " | J_ii = ",    @sprintf("%.2f",p.J_ii),
         " | g_ii = ",    @sprintf("%.2f",p.g_ii),
         " | ν = ",       @sprintf("%.2f",p.scaleg),
         " | I_i = ",     @sprintf("%.2f",p.I_i),
         " | I = ",       @sprintf("%.2f",I_i),
         " | mean(J) = ", @sprintf("%.4f",mean(J)),
         " | std(J) = ",  @sprintf("%.4f",std(J)),
         " | Φ = ",       @sprintf("%s",p.Φname),
         " | λ = ",       @sprintf("%.2f",p.λ),
         " | dist = ",    @sprintf("%s",p.Jdist),
         " | ρ_e = ",     @sprintf("%.2f",p.ρ_e),
         " | dt = ",      @sprintf("%.4f",p.dt),
         " | T = ",       @sprintf("%d",p.T),
         " | s = ",       @sprintf("%d",p.seedJ), "\n")
   printstyled(runparams, color=:green)
end

#################################################
function genfilename(p)
   SUFFIX = string("data/",p.prefix,
                   @sprintf("_N%d",     p.N),
                   @sprintf("_K%d",     p.K),
                   @sprintf("_T%d",     p.T),
                   @sprintf("_d%.4f",   p.dt),
                   @sprintf("_rcDt%.4f",p.rcDt),
                   @sprintf("_tau%.4f", p.τ),
                   @sprintf("_Φ%s",     p.Φname),
                   @sprintf("_λ%.4f",   p.λ),
                   @sprintf("_Λ%s",     p.Jdist),
                   @sprintf("_II%.4f",  p.J_ii),
                   @sprintf("_scJ%.2f", p.scaleJ),
                   @sprintf("_gII%.4f", p.g_ii),
                   @sprintf("_scg%.4f", p.scaleg),
                   @sprintf("_bI%.4f",  p.I_i),
                   @sprintf("_scI%.2f", p.scaleI),
                   @sprintf("_sIC%d",   p.seedIC),
                   @sprintf("_sJ%d",    p.seedJ),
                   ".jld2")
   return SUFFIX
end

#################################################
function writefile(results, p)
   filename = genfilename(p)
   save(filename, "data", results)
end

#################################################
function readfile(p)
   filename = genfilename(p)
   return load(filename)["data"]
end


end
