include("src/model.jl")
using Printf
using JLD2,FileIO

K = round.(Int, 10 .^range(log10(500), log10(20000), length=10))[1:2]

ϕ̄ = zeros(length(K))
f = zeros(length(K))
μ = zeros(length(K))
β = zeros(length(K))

for i = 1:length(K)
   v.dp.K = K[i]
   v.dp.N = K[i]
   v.dp.rcNsamp = min(1000, v.dp.N)

   D,F = v.runsim()
   writefile(D, v.dp)

   p.prefix = "F"
   writefile(F, v.dp)

   ϕ̄[i] = D[:m]
   f[i] = D[:f̄]
   μ[i] = D[:μ̄]
   β[i] = D[:β]
end


