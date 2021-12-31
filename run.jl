include("src/model.jl")
using Printf
using PyPlot
using PyCall
using JLD2,FileIO

K = round.(Int, 10 .^range(log10(500), log10(20000), length=10))[1:7]
scalegs = [-0.25, -0.5]

ϕ̄ = zeros(length(K),2)
f = zeros(length(K),2)
μ = zeros(length(K),2)

for k = 1:2
   v.dp.scaleg = scalegs[k]
   if k==2
      v.dp.J_ii = 1.05
   end
   for i = 1:length(K)
      v.dp.K = K[i]
      v.dp.N = K[i]
      v.dp.rcNsamp = min(1000, v.dp.N)
   
      filename = v.genfilename(v.dp)
   
      if !isfile(filename)
         D,F = v.runsim()
         v.writefile(D, v.dp)
      else
         D = v.readfile(v.dp)
      end
   
      ϕ̄[i,k] = D[:m]
      f[i,k] = D[:f̄]
      μ[i,k] = D[:μ̄]
   end
end
    
lowvar  = "#b2581e"
highvar = "#1e78b2" 

fig, ax = plt.subplots(1,3,figsize=(19/2.54,7/2.54))
mticker = pyimport("matplotlib.ticker")

ax[1].plot(K, ϕ̄[:,1], color=highvar, label="high var")
ax[1].plot(K, ϕ̄[:,2], color=lowvar, linestyle="dashed", label="low var")
ax[1].set(ylabel=L"mean response ($\overline{\phi}$)")
ax[1].legend(frameon=false)

ax[2].plot(K, f[:,1], color=highvar)
ax[2].plot(K, f[:,2], color=lowvar)
ax[2].set(ylabel=L"fraction active ($f$)")

ax[3].plot(K, μ[:,1], color=highvar)
ax[3].plot(K, μ[:,2], color=lowvar)
ax[3].set(ylabel=L"mean active response ($\mu$)")

ax[1].set_xscale("log");ax[1].set_yscale("log")
ax[2].set_xscale("log");ax[2].set_yscale("log")
ax[3].set_xscale("log");ax[3].set_yscale("log")

ax[1].yaxis.set_minor_formatter(mticker.NullFormatter())
ax[2].yaxis.set_minor_formatter(mticker.NullFormatter())
ax[3].yaxis.set_minor_formatter(mticker.NullFormatter())

ax[1].set_yticks([0.01,0.03])
ax[2].set_yticks([0.1,0.3])
ax[3].set_yticks([0.01,0.1])

ax[1].get_yaxis().set_major_formatter(mticker.ScalarFormatter())
ax[2].get_yaxis().set_major_formatter(mticker.ScalarFormatter())
ax[3].get_yaxis().set_major_formatter(mticker.ScalarFormatter())

for i in 1:length(ax)
   ax[i].set_xticks([1000,10000])
   ax[i].set(xlabel=L"in-degree ($K$)")
   ax[i].spines["right"].set_visible(false)
   ax[i].spines["top"].set_visible(false)
   ax[i].xaxis.set_ticks_position("bottom")
   ax[i].yaxis.set_ticks_position("left")
end

tight_layout()
savefig("plots/Kscalings.pdf")


###############################################################
fig, ax = plt.subplots(1,3,figsize=(19/2.54,7/2.54))

v.dp.scaleg = -0.5
v.dp.J_ii   = 1.05
v.dp.N = 1000
v.dp.K = 1000
v.dp.T = 500
D, F = v.runsim()

times = 0 : v.dp.rcDt : (v.dp.T-v.dp.rcDt)
ax[1].plot(times, F[:ϕs][:,1:7], alpha=0.4)
ax[3].hist(v.unroll(F[:xs]), 100, alpha=0.7, density=true, color=lowvar, label="low var")
ax[3].axvline(x=D[:x̄], color=lowvar, linestyle="dashed", alpha=0.7)

v.dp.scaleg = -0.25
v.dp.J_ii   = 2.0
v.dp.N = 1000
v.dp.K = 1000
v.dp.T = 500
D, F = v.runsim()

ax[2].plot(times, F[:ϕs][:,1:7], alpha=0.4)
ax[3].hist(v.unroll(F[:xs]), 100, alpha=0.7, density=true, color=highvar, label="high var")
ax[3].axvline(x=0,color="gray")
ax[3].axvline(x=D[:x̄], color=highvar, linestyle="dashed", alpha=0.7)

ax[1].set(title="low var", xlabel="time",ylabel="response")
ax[2].set(title="high var",xlabel="time",ylabel="response")

ax[1].set_ylim([0,0.6])
ax[2].set_ylim([0,0.6])
ax[3].set_xlim([-4,1])
ax[3].legend(frameon=false)

for i in 1:length(ax)
   ax[i].spines["right"].set_visible(false)
   ax[i].spines["top"].set_visible(false)
   ax[i].xaxis.set_ticks_position("bottom")
   ax[i].yaxis.set_ticks_position("left")
end

tight_layout()
savefig("plots/activities.pdf")

