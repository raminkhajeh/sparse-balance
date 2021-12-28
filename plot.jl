include("src/model.jl")
using PyPlot


fig, ax = plt.subplots(2,3,figsize=(17/2.54,10/2.54))

ax[1].plot(K, ϕ̄)
ax[1].set(ylabel=L"mean response ($\overline{\phi}$)")
ax[1].set(xscale="log",yscale="log")
ax[1].set_yticks([0.01,0.03])

ax[3].plot(K, f)
ax[3].set(ylabel=L"fraction active ($f$)")
ax[3].set(xscale="log",yscale="log")
ax[3].set_yticks([0.1,0.3])

ax[5].plot(K, μ)
ax[5].set(ylabel=L"mean active response ($\mu$)")
ax[5].set(xscale="log",yscale="log")
ax[5].set_yticks([0.02,0.1])

ax[2].plot(K, β)
ax[2].set(ylabel=L"decorrelation rate ($\beta$)")
ax[2].set(xscale="log")

ax[3].hist(,100,alpha=0.5)
ax[3].hist(,100,alpha=0.5)

for i = 1:length(ax)
   ax[i].set_xticks([1000,10000])
   ax[i].set(xlabel=L"in-degree ($K$)")
   ax[i].spines["right"].set_visible(false)
   ax[i].spines["top"].set_visible(false)
   ax[i].xaxis.set_ticks_position("bottom")
   ax[i].yaxis.set_ticks_position("left")
end

tight_layout()
savefig("plots/fig.pdf")
