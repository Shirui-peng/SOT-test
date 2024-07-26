include("../src/SOT.jl")
using .SOT, Printf, LinearAlgebra, Dates
using HDF5, DataFrames, CSV, Statistics, StatsBase
using Random, Distributions, SparseArrays
using PyPlot, Interpolations, Cubature, LsqFit
PyPlot.matplotlib[:rc]("mathtext",fontset="cm")        #computer modern font 
PyPlot.matplotlib[:rc]("font",family="STIXGeneral")
rc("font", size=8)
rc("axes", titlesize="medium")

### T waves
h11a0 = [141.96,38.58]
#te = [DateTime(2018, 11, 24)]
te = DateTime(2000, 1, 1) .+ Millisecond.(h5read("data/japan_H11N3_x2y_16ex.h5", "t"))
m = length(te)
#lone,late,azme = ones(m)*h11a0[1],ones(m)*h11a0[2],zeros(m)
lone = h5read("data/japan_H11N3_x2y_16ex.h5", "lon")
late = h5read("data/japan_H11N3_x2y_16ex.h5", "lat")
azme = h5read("data/japan_H11N3_x2y_16ex.h5", "θ")
nms = ["longitude","latitude","time","azimuth"]
events = DataFrame([lone late te azme], nms)

station = "H11N3"
nexclude = 16
ppfile = @sprintf("../results/pairs/japan_%s_ppairs_2.5a3.5hz_%dex.csv",station,nexclude)
tpvfile = @sprintf("../results/pairs/japan_%s_tpairs_2.5a3.5hz_%dex.csv",station,nexclude)
tstation = [166.90986,19.71786]

# manually exclude pairs
tpairs = CSV.read(tpvfile, DataFrame)
ppairs = CSV.read(ppfile, DataFrame)

#tpairs.Δτt = tpairs.Δτ .- tpairs.Δτp 
# collect delays into data vector
ppairs = innerjoin(ppairs, select(tpairs, [:event1, :event2]), on=[:event1, :event2])
ppairs.stn6 = [s[1:6] for s in ppairs.station]
unique!(ppairs, [:stn6,:event1,:event2])
select!(ppairs, Not([:stn6,:cc,:n1,:n2]))

evtpos = [142.85,38.10]

pstations = ["IU.MAJO.00.BHZ","IU.MAJO.10.BHZ","PS.TSK..BHZ", "II.ERM.00.BHZ", "G.INU.00.BHZ"]
pstnlats,pstnlons = [36.55,36.55,36.21,42.02,35.35],[138.2,138.2,140.11,143.16,137.03]
pstations = DataFrame(station=pstations,slat=pstnlats,slon=pstnlons)

nt1 = size(tpairs, 1)
np1 = size(ppairs, 1)

σtrend = 0.01

# annual cycle prior (s)
σannual = 0.1

# semi-annual cycle prior (s)
σsemiannual = 0.1

λt,στ,σp,σx,σn,σnp,σh = 62,0.29,0.97,0.019,8.1e-3,2.2e-3,0.029
x = [λt,στ,σp,σx,σn,σnp,σh]

θgrid = -1:0.1:20
cgrid = SOT.cint.(θgrid)
itp_c = linear_interpolation(θgrid, cgrid, extrapolation_bc=0)
cgrid = SOT.tcint.(θgrid)
itp_tc = linear_interpolation(θgrid, cgrid, extrapolation_bc=0)
3.7σtrend/sqrt(itp_tc(0))

θgrid = -10:0.2:10
xgrid,ygrid = 135:0.1:170,15:0.1:45
ctgrid = SOT.cinta.(xgrid,ygrid',reshape(θgrid,1,1,:))
nodes = (xgrid, ygrid, θgrid)
itp_ca = interpolate(nodes, ctgrid, Gridded(Linear()))
tctgrid = SOT.tcinta.(xgrid,ygrid',reshape(θgrid,1,1,:))
nodes = (xgrid, ygrid, θgrid)
itp_tca = interpolate(nodes, tctgrid, Gridded(Linear()))
println("covariance integration finished.")

c2grid = SOT.cint2da.(xgrid,ygrid')
tc2grid = SOT.tcint2da.(xgrid,ygrid')
shape = "rec"
rc2grid = SOT.cint2da.(xgrid,ygrid';shape)
rtc2grid = SOT.tcint2da.(xgrid,ygrid';shape)

rc("font", size=8)
rc("axes", titlesize="medium")
fig,axs = plt.subplots(2,2,figsize=(6.4,4.8),sharex=true,sharey=true)
axs = reshape(axs,(4,))
axs[1].pcolormesh(xgrid,ygrid,ctgrid[:,:,argmin(abs.(θgrid))]',cmap="Blues",rasterized=true)
#axs[2].pcolormesh(xgrid,ygrid,c2grid',cmap="Blues",rasterized=true)
axs[3].pcolormesh(xgrid,ygrid,tctgrid[:,:,argmin(abs.(θgrid))]',cmap="Blues",rasterized=true)
#axs[5].pcolormesh(xgrid,ygrid,tc2grid',cmap="Blues",rasterized=true)
axs[2].pcolormesh(xgrid,ygrid,rc2grid',cmap="Blues",rasterized=true)
axs[4].pcolormesh(xgrid,ygrid,rtc2grid',cmap="Blues",rasterized=true)
axs[1].set_ylabel("latitude")
axs[2].set_ylabel("latitude")
#axs[3].set_ylabel("latitude")
axs[2].set_xlabel("longitude")
axs[4].set_xlabel("longitude")
axs[1].set_title("stochastic covariance")
axs[3].set_title("trend covariance")
fig.tight_layout()
fig.savefig("results/cgrid.pdf")

tstation = [166.652,19.283]
θ1,θ2 = -6.5,3.5
ctdfile = "results/ctds_kuroshio.csv"
ctds = CSV.read(ctdfile, DataFrame)
θ0 = SOT.azimuth(tstation[2],tstation[1],evtpos[2],evtpos[1])
idxa = θ1 .<= SOT.azimuth.(tstation[2],tstation[1],ctds.y,ctds.x) .- θ0 .<= θ2
ctds = ctds[vec(idxa),:]
ctds = ctds[(ctds.y .<=36) .| (ctds.x .>= 141),:]
ctds = ctds[(DateTime(2001,2,25) .≤ ctds.t .< DateTime(2021,8,1)),:]

argofile = "results/japan_1997-01_2022-12.csv"

# manually exclude pairs
floats = CSV.read(argofile, DataFrame)

t0,t1,y0,z0 = DateTime(2001,2,1),DateTime(2021,8,1),30,-1.9e3
floats = floats[(floats.z0.<z0),:]
floats = floats[(t0 .≤ floats.t .< t1),:]
floats = floats[.!(isnan.(floats.Δsp1)),:]
floats = floats[(floats.y .<=36) .| (floats.x .>= 141),:]
unique!(floats, [:x,:y,:t])
#floats = floats[(floats.y .<=34) .| (floats.y .>= 36),:]
#floats = floats[(floats.y .>=30) .& (floats.x .>= 140),:]

d0 = 300
Δt = 120
nxk = 321
#nearfloats,ictds = SOT.nearctds(floats,ctds,d0,Δt)
idxa = θ1.<= SOT.azimuth.(tstation[2],tstation[1],floats.y, floats.x) .- θ0 .<= θ2
nearfloats = floats[vec(idxa),:]
#ctds = ctds[Array{Int}(ictds),:]
println(size(nearfloats,1))
println(size(ctds,1))
#7555 1223

xtrd=3.7
Pa,xa,Pc,xc,Pac,xac,C0,dC = SOT.invertctd(xtrd*σtrend, σannual,ctds,nearfloats,100/sqrt(2),60,0.44,0.44,itp_tc)
m = size(ctds,1)
trc = Dates.value.(ctds.t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
tm = (max(trc...)+min(trc...))/2
m = length(trc)
ω = 2π/SOT.meanyear
D = [I(m) I(m) Diagonal(trc.-tm) cos.(ω*trc) sin.(ω*trc) cos.(2ω*trc) sin.(2ω*trc)]

tbnd = xtrd*σtrend/sqrt(itp_tc(0))/SOT.meanyear
dxc = SOT.dist.(ctds.x', ctds.y', ctds.x, ctds.y)/1e3
U,S,Vt = svd(Symmetric(exp.(-dxc/200)))
zit = Diagonal(sqrt.(S).^-1)*U'*(xc[1+m:2m])/(0.47/2)
println(sum(abs.(zit).<2)/m)
zit = Diagonal(sqrt.(S).^-1)*U'*(xc[1+2m:3m])/tbnd
println(sum(abs.(zit).<2)/m)

dP = 2*(Pa+Pc-C0)+dC+dC'
@printf("largest entry: %.2e\n",max(abs.(dP-2*Pac)...))
CiP = cholesky(Symmetric(inv(D*(Pa+Pc-dP)*D')))
ziP = CiP.U*D*(xa.-xc)
#U,S,Vt = svd(Symmetric(Pa+Pc-2Pac))
#ziP = Diagonal(sqrt.(S).^-1)*U'*(xa.-xc)
frac = sum(-2 .< ziP .< 2)/length(D*xa)
println(frac) #86%

p = 0.5:0.5:99.5
α = 0.05
zτh = h5read("results/argotris_H11.h5", "zτ")
zτw = h5read("results/argotris_WAKE.h5", "zτ")
ziP = h5read("results/argotris_ctd.h5", "zs")
println([std(zτh),std(zτw),std(ziP)])
qsc = percentile(ziP,p)
qs4,qs5 = percentile(zτh,p),percentile(zτw,p)
qt = quantile.(Normal(), p/100)
lb,ub = zeros(length(p),3),zeros(length(p),3)
for i = 1:length(p)
    for (j,nt) in enumerate([length(ziP),length(zτh),length(zτw)])
        d = NoncentralT(nt-1,-sqrt(nt)*qt[i])
        lb[i,j],ub[i,j] = -quantile(d, 1-α/2)/sqrt(nt),-quantile(d, α/2)/sqrt(nt)
    end
end
rc("font", size=10)
fig,ax=subplots(3,1,figsize=(4.8,4.8),sharex=true)
c3 = ["#1b9e77","#d95f02","#7570b3"]
ax[1].set_title("H11")
ax[1].set_title("(a)",loc="left")
ax[1].plot(qt,qs4-qt,label="H11",c=c3[1])
ax[1].fill_between(qt, (lb[:,2]-qt), (ub[:,2]-qt), alpha=.2, zorder=3, color=c3[1], linewidth=0)
ax[2].set_title("WAKE")
ax[2].set_title("(b)",loc="left")
ax[2].plot(qt,qs5-qt,label="WAKE",c=c3[2])
ax[2].fill_between(qt, (lb[:,3]-qt), (ub[:,3]-qt), alpha=.2, zorder=3, color=c3[2], linewidth=0)
ax[3].set_title("Shipboard CTD")
ax[3].set_title("(c)",loc="left")
ax[3].plot(qt,qsc-qt,label="shipboard CTD",c=c3[3])
ax[3].fill_between(qt, (lb[:,1]-qt), (ub[:,1]-qt), alpha=.2, zorder=3, color=c3[3], linewidth=0)
#ax.axhline(0,color="black",ls=":",lw=1)
ax[3].set_xlabel("theoretical quantile")
ax[2].set_ylabel("sample quantile \$-\$ theoretical quantile")
#fig.text(0.01, 0.5, "sample quantile \$-\$ theoretical quantile", va="center", rotation="vertical")
#ax.legend(frameon=false)
ax[3].set_xlim([qt[1],qt[end]])
fig.tight_layout()
fig.savefig("results/quantile3_kuroshio.pdf",bbox_inches="tight",dpi=300)
println(sum(lb[:,2].<=qs4.<=ub[:,2])/length(qt))
println(sum(lb[:,3].<=qs5.<=ub[:,3])/length(qt))
println(sum(lb[:,1].<=qsc.<=ub[:,1])/length(qt))

h5open("results/argotris_ctd.h5", "w") do file
    write(file, "sa", D*xa)
    write(file, "sc", D*xc)
    write(file, "sac", D*xac)
    write(file, "zs", ziP)
end

nxk = 321
is5,in5 = argmin(azme),argmin(azme.-5)
xks, yks = SOT.findpath(tstation, (lone[is5],late[is5]), nxk)
#xkn, ykn = SOT.findpath(tstation, (lone[in5],late[in5]), nxk)
d0 = 300
idxs = minimum(SOT.dist.(xks', yks', floats.x, floats.y),dims=2)./1e3 .< d0
#idxn = minimum(SOT.dist.(xkn', ykn', floats.x, floats.y),dims=2)./1e3 .< d0
θ0 = SOT.azimuth(tstation[2],tstation[1],evtpos[2],evtpos[1])
idxa = min(azme...) .<= SOT.azimuth.(tstation[2],tstation[1],floats.y, floats.x) .- θ0 .<= max(azme...)
floats = floats[vec(idxa) .| vec(idxs),:]
println(size(floats,1))

tm = Dates.value.(t0+(t1-t0)/2 - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
trd = Dates.value.(floats.t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
X = [ones(size(floats,1)) (trd.-tm)]
β = (X'*X)\(X'*floats.Δsp1)

ω = 2π/SOT.meanyear
#Ξ = Diagonal(([1.5σtrend/SOT.meanyear; σannual*ones(4)/sqrt(2)]/3.2).^2)
iN = (σs[1]^2+σn[1]^2)^-1*I
lata,lona = 20:40,140:165
npa = zeros(length(lata),length(lona))
ma,trenda = zeros(length(lata),length(lona)),zeros(length(lata),length(lona))
for (i,lat) in enumerate(lata)
  for (j,lon) in enumerate(lona)
    idx = (lat-3 .≤ floats.y .< lat+3) .& (lon-3 .≤ floats.x .< lon+3)
    profiles = floats[idx,:]
    npa[i,j] = size(profiles,1)
    if npa[i,j]>10
        println(npa[i,j])
        trd = Dates.value.(profiles.t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
        X = [ones(size(profiles,1)) (trd.-tm)]
        β = (X'*X)\(X'*profiles.Δsp1)
        trenda[i,j] = β[2]*SOT.meanyear*3.2
        ma[i,j] = β[1]
    end
  end
end

fig,ax=subplots()
im = ax.pcolormesh(lona,lata,ma,cmap="RdBu_r", vmin=-.5, vmax=.5, shading="nearest", rasterized="true")
fig.colorbar(im,ax=ax)
fig.tight_layout()
fig.savefig("results/kuroshio_argo_mean.pdf")

fig,axs=subplots(1,3,figsize=(5.4,4.8),sharey=true)
ax=axs[1]
ax.plot(1e3trenda,lata)
ax.fill_betweenx(lata, 1e3(trenda-2ea), 1e3(trenda+2ea), alpha=.2, zorder=3, color="tab:blue", linewidth=0)
ax.set_xlabel("trend (ms/yr)")
ax.set_ylabel("latitude")
axs[2].plot(npa,lata)
axs[2].set_xlabel("profile number")
cor = autocor(trenda)
xdata = range(0, stop=length(cor)-1, length=length(cor))
fit = curve_fit(model, xdata, cor, [2.0])
cor0 = [cor[2:14][end:-1:1]; cor[1:14]]
axs[3].plot(cor0,lata,label="data")
axs[3].plot(exp.(-abs.(lata.-30)./coef(fit)[1]),lata,label=@sprintf("exp fit, %.1f\$^\\circ\$",coef(fit)[1]))
axs[3].legend(frameon=false)
axs[3].set_xlabel("correlation")
fig.tight_layout()
fig.savefig("results/kuroshio_argo_lattrend.pdf")

lona = 136:170
npa,trenda,ea = zeros(length(lona)),zeros(length(lona)),zeros(length(lona))
for (i,lon) in enumerate(lona)
  profiles = floats[(lon-.5 .≤ floats.x .< lon+.5),:]
  npa[i] = size(profiles,1)
  trd = Dates.value.(profiles.t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
  X = [(trd.-tm) cos.(ω*trd) sin.(ω*trd) cos.(2ω*trd) sin.(2ω*trd)]
  P = inv(inv(Matrix(Ξ))+X'*iN*X)
  β = P*(X'*iN*profiles.Δsp1)
  trenda[i] = β[1]*SOT.meanyear*3.2
  ea[i] = sqrt(P[1])*SOT.meanyear*3.2
end

fig,axs=plt.subplots(3,1,figsize=(5.4,4.8),sharex=true)
ax=axs[1]
ax.plot(lona,1e3trenda)
ax.fill_between(lona, 1e3(trenda-2ea), 1e3(trenda+2ea), alpha=.2, zorder=3, color="tab:blue", linewidth=0)
ax.set_ylabel("trend (ms/yr)")
axs[3].set_xlabel("longitude")
axs[2].plot(lona,npa)
axs[2].set_ylabel("profile number")
cor = autocor(trenda)
xdata = range(0, stop=length(cor)-1, length=length(cor))
@. model(x, p) = exp(-x/p[1])
fit = curve_fit(model, xdata, cor, [2.0])
cor0 = [cor[end:-1:2]; cor]
axs[3].plot(lona[3:end-2],cor0,label="data")
axs[3].plot(lona,exp.(-abs.(lona.-153)./coef(fit)[1]),label=@sprintf("exp fit, %.1f\$^\\circ\$",coef(fit)[1]))
axs[3].legend(frameon=false)
axs[3].set_ylabel("correlation")
fig.tight_layout()
fig.savefig("results/kuroshio_argo_lontrend.pdf")

λt = [62,57]

λx = [65,59]

# scale (s/Mm)
σs = [0.49,0.32]

# noise (s/Mm)
σn = [0.041,0.016]

tstation = [166.652,19.283]
xrg = h5read("results/kuroshio_rg.h5", "x")
yrg = h5read("results/kuroshio_rg.h5", "y")
nxrg,nyrg=length(xrg),length(yrg)
a0 = SOT.azimuth.(tstation[2],tstation[1],evtpos[2],evtpos[1])
arg = SOT.azimuth.(tstation[2],tstation[1],yrg',xrg) .- a0
drg = SOT.dist.(tstation[1], tstation[2], xrg, yrg')/1e6
aidx = abs.(drg.*sind.(arg)) .<= 3.2*tand(5)
Δsrg = h5read("results/kuroshio_rg.h5", "Δsp1")[:,:,1]
aidx = aidx .& .!iszero.(Δsrg) .& (drg.<=3.5)
h5fl = "results/kuroshio_mt_H11.h5"
xa,xt,xat = h5read(h5fl, "xa"),h5read(h5fl, "xt"),h5read(h5fl, "xat")
Pa,Pt,Pat = h5read(h5fl, "Pa"),h5read(h5fl, "Pt"),h5read(h5fl, "Pat")
gmn,gtd = fill(NaN, nxrg, nyrg, 2,2),fill(NaN, nxrg,nyrg,2,3)
ng=1
for i = 1:nxrg, j = 1:nyrg
    if aidx[i,j]
        gmn[i,j,1,1],gmn[i,j,1,2]=xa[ng],xat[ng]
        gtd[i,j,1,1],gtd[i,j,1,2],gtd[i,j,1,3] = xa[sum(aidx)+ng],xt[sum(aidx)+ng],xat[sum(aidx)+ng]
        gmn[i,j,2,1],gmn[i,j,2,1]=sqrt(diag(Pa)[ng]),sqrt(diag(Pat)[ng])
        gtd[i,j,2,1],gtd[i,j,2,2],gtd[i,j,2,3]=sqrt(diag(Pa)[sum(aidx)+ng]),sqrt(diag(Pt)[sum(aidx)+ng]),sqrt(diag(Pat)[sum(aidx)+ng])
        global ng+=1
    end
end

fig,ax=subplots(1,2,figsize=(190/25.4, 190/25.4/3),sharex=true,sharey=true)
im1=ax[1].pcolormesh(xrg,yrg,gmn[:,:,1,1]',vmin=-1,vmax=1,cmap="RdBu",shading="auto",rasterized=true)
ax[1].set_title("bias mean")
ax[1].set_xlabel("longitude")
ax[1].set_ylabel("latitude")
ax[2].set_title("bias uncertainty")
im2=ax[2].pcolormesh(xrg,yrg,vmin=0,vmax=0.5,gmn[:,:,2,1]',shading="auto",rasterized=true)
ax[2].set_xlabel("longitude")
fig.tight_layout()
cbaxes = fig.add_axes([0.09, 0.28, 0.05, 0.02]) 
cbaxes.text(1.2, -3.5,"(s/Mm)")
cbar = fig.colorbar(im1, cax=cbaxes, ticks=[-1,0,1], orientation="horizontal")
cbaxes = fig.add_axes([0.56, 0.28, 0.05, 0.02]) 
cbaxes.text(0.65, -3.5,"(s/Mm)")
cbar = fig.colorbar(im2, cax=cbaxes, ticks=[0,0.5], orientation="horizontal")
fig.savefig("results/kuroshio_mean_H11.pdf")

fig,ax=subplots(3,2,figsize=(190/25.4, 190/25.4),sharex=true,sharey=true)
im1=ax[1,1].pcolormesh(xrg,yrg,1e3SOT.meanyear*gtd[:,:,1,1]',vmin=-40,vmax=40,
                      cmap="RdBu",shading="auto",rasterized=true)
ax[1,1].set_title("point trend mean")
ax[1,2].set_title("point trend uncertainty")
im=ax[2,1].pcolormesh(xrg,yrg,1e3SOT.meanyear*gtd[:,:,1,2]',vmin=-40,vmax=40,
                      cmap="RdBu",shading="auto",rasterized=true)
ax[2,1].set_title("SOT trend mean")
ax[2,2].set_title("SOT trend uncertainty")
im=ax[3,1].pcolormesh(xrg,yrg,1e3SOT.meanyear*gtd[:,:,1,3]',vmin=-40,vmax=40,
                      cmap="RdBu",shading="auto",rasterized=true)
ax[3,1].set_title("full trend mean")
ax[3,2].set_title("full trend uncertainty")
im2=ax[1,2].pcolormesh(xrg,yrg,1e3SOT.meanyear*gtd[:,:,2,1]',vmin=13,vmax=20,
                       shading="auto",rasterized=true)
im=ax[2,2].pcolormesh(xrg,yrg,1e3SOT.meanyear*gtd[:,:,2,2]',vmin=13,vmax=20,
                      shading="auto",rasterized=true)
im=ax[3,2].pcolormesh(xrg,yrg,1e3SOT.meanyear*gtd[:,:,2,3]',vmin=13,vmax=20,
                      shading="auto",rasterized=true)
for i=1:2,j=1:3
    k = 2*(j-1)+i
    ax[j,i].set_title("($(('a':'z')[k]))", loc="left")
    ax[3,i].set_xlabel("longitude")
    ax[j,1].set_ylabel("latitude")
end
fig.tight_layout()
cbaxes = fig.add_axes([0.09, 0.72, 0.05, 0.01]) 
cbaxes.text(58, -2.3,"(ms/Mm/yr)")
cbar = fig.colorbar(im1, cax=cbaxes, ticks=[-40,0,40], orientation="horizontal")
cbaxes = fig.add_axes([0.56, 0.72, 0.05, 0.01]) 
cbaxes.text(22, -2.3,"(ms/Mm/yr)")
cbar = fig.colorbar(im2, cax=cbaxes, ticks=[14,17,20], orientation="horizontal")
fig.savefig("results/kuroshio_trend_H11.pdf")

gxy = zeros(sum(aidx),2)
na = 1
for i = 1:nxrg, j = 1:nyrg
    if aidx[i,j] 
        gxy[na,:]=[xrg[i],yrg[j]]
        na+=1
    end
end

floats = vcat(CSV.read(ctdfile, DataFrame)[:,1:6],CSV.read(argofile, DataFrame))
floats = floats[(floats.y .<=36) .| (floats.x .>= 141),:]
floats = floats[.!(isnan.(floats.Δsp1)),:]
unique!(floats, [:x,:y,:t])
nearfloats = SOT.nearargo(floats,events,nxk,d0,Δt,tstation)

xrg = h5read("results/kuroshio_rg.h5", "x")
yrg = h5read("results/kuroshio_rg.h5", "y")
xrg,yrg = xrg[1]:0.5:xrg[end],yrg[1]:0.5:yrg[end]
nxrg,nyrg=length(xrg),length(yrg)
a0 = SOT.azimuth.(tstation[2],tstation[1],evtpos[2],evtpos[1])
arg = SOT.azimuth.(tstation[2],tstation[1],yrg',xrg) .- a0
drg = SOT.dist.(tstation[1], tstation[2], xrg, yrg')/1e6
aidx = abs.(drg.*sind.(arg)) .<= 3.2*tand(5)
aidx = aidx .& (drg.<=3.2)
gxy = zeros(sum(aidx),3)
ted = Dates.value.(te[1] - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
na = 1
for i = 1:nxrg, j = 1:nyrg
    if aidx[i,j] 
        gxy[na,:]=[xrg[i],yrg[j],ted]
        na+=1
    end
end

Pa,xa,Pt,xt,Pat,xat = SOT.invertargotag(gxy,x,σtrend, σannual, σsemiannual,tpairs,ppairs,tstation,evtpos,pstations,nearfloats,λx[1],σs[1],σn[1],itp_c,itp_ca,itp_tc,itp_tca;tλx=200)

gae = fill(NaN, nxrg, nyrg, 3,2)
ng=1
for i = 1:nxrg, j = 1:nyrg
    if aidx[i,j]
        gae[i,j,1,1],gae[i,j,1,2]=xa[ng],sqrt(diag(Pa)[ng])
        gae[i,j,2,1],gae[i,j,2,2]=xt[ng],sqrt(diag(Pt)[ng])
        gae[i,j,3,1],gae[i,j,3,2]=xat[ng],sqrt(diag(Pat)[ng])
        ng+=1
    end
end
xk0, yk0 = SOT.findpath(tstation, h11a0, nxk)
nodes = (xrg, yrg)
fig,ax=subplots(2,1,figsize=(6.4,4.8),sharex=true)
for i=1:3
    for j=1:2
        itp = interpolate(nodes, gae[:,:,i,j], Gridded(Linear()))
        ak = itp.(xk0,yk0)
        println(size(gae[:,:,i,j]))
        ax[j].plot(0:10:3200,ak)
    end
end
fig.savefig("results/kuroshio_anmly0d_H11.pdf")

fig,ax=subplots(3,2,figsize=(190/25.4, 190/25.4),sharex=true,sharey=true)
im1=ax[1,1].pcolormesh(xrg,yrg,gae[:,:,1,1]',vmin=-1,vmax=1,
                       cmap="RdBu",shading="auto",rasterized=true)
ax[1,1].plot(xk0,yk0,lw=1,c="k")
ax[1,2].plot(xk0,yk0,lw=1,c="k")
for i=2:3
    ax[i,1].pcolormesh(xrg,yrg,gae[:,:,i,1]',vmin=-1,vmax=1,
                       cmap="RdBu",shading="auto",rasterized=true)
    ax[i,1].plot(xk0,yk0,lw=1,c="k")
    ax[i,2].plot(xk0,yk0,lw=1,c="k")
end
ax[1,1].set_title("point anomaly mean")
ax[2,1].set_title("\$T\$-wave anomaly mean")
ax[2,1].set_title("full anomaly mean")
vem = 0.3
im2=ax[1,2].pcolormesh(xrg,yrg,gae[:,:,1,2]',vmin=vem,vmax=0.5,
                       shading="auto",rasterized=true)
for i=2:3
    ax[i,2].pcolormesh(xrg,yrg,gae[:,:,i,2]',vmin=vem,vmax=0.5,
                       shading="auto",rasterized=true)
end
ax[1,2].set_title("point anomaly uncertainty")
ax[2,2].set_title("\$T\$-wave anomaly uncertainty")
ax[3,2].set_title("full anomaly uncertainty")
for i=1:2,j=1:3
    k = 2*(j-1)+i
    ax[j,i].set_title("($(('a':'z')[k]))", loc="left")
    ax[3,i].set_xlabel("longitude")
    ax[j,1].set_ylabel("latitude")
end
fig.tight_layout()
cbaxes = fig.add_axes([0.09, 0.72, 0.05, 0.01]) 
cbaxes.text(1.3, -2.3,"(s/Mm)")
cbar = fig.colorbar(im1, cax=cbaxes, ticks=[-1,0,1], orientation="horizontal")
cbaxes = fig.add_axes([0.56, 0.72, 0.05, 0.01]) 
cbaxes.text(0.57, -2.3,"(s/Mm)")
cbar = fig.colorbar(im2, cax=cbaxes, ticks=[vem,0.5], orientation="horizontal")
fig.savefig("results/kuroshio_anomaly_H11.pdf")

x0 = log.([66,65,0.5,0.05])
f(x) = -1 .* SOT.loglikelihoodargo(x, floats.Δsp1-X*β, floats.t, floats.x, floats.y;dyear=2)
#f(x) = -1 .* SOT.loglikelihoodargo(x, y, [nearfloats.t;ctds.t], [nearfloats.x;ctds.x], [nearfloats.y;ctds.y];σtrend, σannual, σsemiannual,itp_tc,full=true)
xk,Hk = SOT.BFGS(f,x0,1e-5,100; xrtol=1e-6)
@printf("MLE: %s\n",exp.(xk))
@printf("lb: %s\n ub: %s\n",exp.(xk.-2*sqrt.(diag(Hk))),exp.(xk.+2*sqrt.(diag(Hk))))
MLE: [66.43355053937229, 65.01079879224713, 0.4433527693465209, 0.03278997159787347]
lb: [59.073668614374654, 62.10862382500905, 0.4242714809112686, 0.023603940395737964]
ub: [74.71038689128233, 68.04858487146532, 0.4632922242735847, 0.04555096392225631]

### 0 deg 
Pa,xa,Pt,xt,Pat,xat = SOT.invertargot0(te,azme,x,2.5σtrend, σannual, σsemiannual,tpairs,ppairs,tstation,evtpos,pstations,nearfloats,λx[1],σs[1],σn[1],itp_c,itp_ca,itp_tc,itp_tca)
#Pa,xa,Pt,xt,Pat,xat = SOT.invertargot0(te,azme,x,2.5σtrend, σannual, σsemiannual,tpairs,ppairs,tstation,evtpos,pstations,nearfloats,λx[1],σs[1],σn[1],itp_c,itp_ca)

ter = Dates.value.(te - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
h5name = "results/argo_H11_0d.h5"
xa = h5read(h5name, "xa")
xt = h5read(h5name, "xt")
xat = h5read(h5name, "xat")
Pa = h5read(h5name, "Pa")
Pt = h5read(h5name, "Pt")
Pat = h5read(h5name, "Pat")
ter = h5read(h5name, "t")
te = DateTime(2000, 1, 1) .+ Millisecond.(ter)
m = length(te)
ter /= (1000*24*3600)
tm = (max(ter...)+min(ter...))/2
ω = 2π/SOT.meanyear
D = [I(m) Diagonal(ter.-tm) cos.(ω*ter) sin.(ω*ter) cos.(2ω*ter) sin.(2ω*ter)]
τa,τt,τat = D*xa,D*xt,D*xat
ea,et,eat = sqrt.(diag(D*Pa*D')),sqrt.(diag(D*Pt*D')),sqrt.(diag(D*Pat*D'))
println(1-sqrt(sum(eat.^2)/sum(ea.^2)))
println(1-sqrt(diag(Pat)[m+1]/diag(Pa)[m+1]))
@printf("Argo: trend %.1e s/yr, etrend %.1e s/yr\n",xa[m+1]*SOT.meanyear,2sqrt(diag(Pa)[m+1])*SOT.meanyear)
@printf("Twave: trend %.1e s/yr, etrend %.1e s/yr\n",xt[m+1]*SOT.meanyear,2sqrt(diag(Pt)[m+1])*SOT.meanyear)
@printf("joint: trend %.1e s/yr, etrend %.1e s/yr\n",xat[m+1]*SOT.meanyear,2sqrt(diag(Pat)[m+1])*SOT.meanyear)
K = -6
@printf("Argo: trend %.1e K/yr, etrend %.1e K/yr\n",xa[m+1]*SOT.meanyear/K,-2sqrt(diag(Pa)[m+1])*SOT.meanyear/K)
@printf("Twave: trend %.1e K/yr, etrend %.1e K/yr\n",xt[m+1]*SOT.meanyear/K,-2sqrt(diag(Pt)[m+1])*SOT.meanyear/K)
@printf("joint: trend %.1e K/yr, etrend %.1e K/yr\n",xat[m+1]*SOT.meanyear/K,-2sqrt(diag(Pat)[m+1])*SOT.meanyear/K)
# 44(15)
#Argo: trend -4.3e-02 s/yr, etrend 2.7e-02 s/yr
#Twave: trend -3.9e-02 s/yr, etrend 2.5e-02 s/yr
#joint: trend -4.3e-02 s/yr, etrend 2.3e-02 s/yr

nxk = 321
xk,yk = SOT.findpath(h11a0, tstation, nxk)
xrg = h5read("results/kuroshio_rg.h5", "x")
yrg = h5read("results/kuroshio_rg.h5", "y")
trg = DateTime(2000, 1, 1) .+ Millisecond.(h5read("results/kuroshio_rg.h5", "t"))
Δsrg = h5read("results/kuroshio_rg.h5", "Δsp1")
nodes = (xrg, yrg)
ntrg = length(trg)
τrg = zeros(ntrg)
for i = 1:ntrg
    itprg = interpolate(nodes, Δsrg[:,:,i], Gridded(Linear()))
    τrg[i] = sum(itprg.(xk,yk))*3.2/nxk
end

K = -6
rc("font", size=8)
rc("axes", titlesize="medium")
fig,axs=subplots(2,1,figsize=(190/25.4,190/25.4/2),sharex=true)
ax = axs[1]
p1, = ax.plot(te,τa/K, linewidth=1,color="tab:blue")
ax.plot(trg,τrg/K, linewidth=1, linestyle="--",color="tab:blue")
ax.fill_between(te, (τa-2ea)/K, (τa+2ea)/K, alpha=.2, zorder=3, color="tab:blue", linewidth=0)
p2, = ax.plot(te,τt/K, linewidth=1,color="tab:orange")
ax.fill_between(te, (τt-2et)/K, (τt+2et)/K, alpha=.2, zorder=3, color="tab:orange", linewidth=0)
p3, = ax.plot(te,τat/K, linewidth=1,color="tab:green")
ax.fill_between(te, (τat-2eat)/K, (τat+2eat)/K, alpha=.2, zorder=3, color="tab:green", linewidth=0)
ax.legend([p1, p2, p3], ["point mode 1 ","\$T\$ waves", "point with \$T\$ waves"], ncol=3, loc="upper center", frameon=false)
ax.set_ylabel("temperature anomaly (K)")
ax.set_title("\$43\\%\$(\$13\\%\$) full(trend) RMSE reduction for point with \$T\$ waves, \$\\alpha={0^\\circ}\$")
ax = axs[2]
ax.plot(te,2ea/abs(K),color="tab:blue", linewidth=1)
ax.plot(te,2et/abs(K),color="tab:orange", linewidth=1)
ax.plot(te,2eat/abs(K),color="tab:green", linewidth=1)
ax.set_ylabel("error (K)")
ax.set_xlim([DateTime(2007,7,1),DateTime(2022,2,1)])
fig.align_ylabels()
fig.tight_layout()
fig.savefig("results/H11_argo_0deg_tc.pdf")

fig,ax=subplots()
ax.plot(te,τa-τt)
ax.set_xlim([DateTime(2018,1,1),DateTime(2020,7,1)])
fig.tight_layout()
fig.savefig("results/H11_argo_0deg_diff.pdf")

#te[argmax(τa-τt)] 2018-11-24

### subbasin
h5name = "results/argo_H11_2d.h5"
xa = h5read(h5name, "xa")
xt = h5read(h5name, "xt")
xat = h5read(h5name, "xat")
Pa = h5read(h5name, "Pa")
Pt = h5read(h5name, "Pt")
Pat = h5read(h5name, "Pat")
ter = h5read(h5name, "t")
te = DateTime(2000, 1, 1) .+ Millisecond.(ter)
ter /= (1000*24*3600)
m = length(te)
tm = (max(ter...)+min(ter...))/2
ω = 2π/SOT.meanyear
D = [I(m) (ter.-tm) cos.(ω*ter) sin.(ω*ter) cos.(2ω*ter) sin.(2ω*ter)]
τa,τt,τat = D*xa,D*xt,D*xat
ea,et,eat = sqrt.(diag(D*Pa*D')),sqrt.(diag(D*Pt*D')),sqrt.(diag(D*Pat*D'))
println(1-sqrt(sum(eat.^2)/sum(ea.^2)))
println(1-sqrt(diag(Pat)[m+1]/diag(Pa)[m+1]))
K = -6
@printf("Argo: trend %.1e K/yr, etrend %.1e K/yr\n",xa[m+1]*SOT.meanyear/K,2sqrt(diag(Pa)[m+1])*SOT.meanyear/K)
@printf("Twave: trend %.1e K/yr, etrend %.1e K/yr\n",xt[m+1]*SOT.meanyear/K,2sqrt(diag(Pt)[m+1])*SOT.meanyear/K)
@printf("joint: trend %.1e K/yr, etrend %.1e K/yr\n",xat[m+1]*SOT.meanyear/K,2sqrt(diag(Pat)[m+1])*SOT.meanyear/K)

xrg = h5read("results/kuroshio_rg.h5", "x")
yrg = h5read("results/kuroshio_rg.h5", "y")
trg = DateTime(2000, 1, 1) .+ Millisecond.(h5read("results/kuroshio_rg.h5", "t"))
Δsrg = h5read("results/kuroshio_rg.h5", "Δsp1")

a0 = SOT.azimuth.(tstation[2],tstation[1],evtpos[2],evtpos[1])
arg = SOT.azimuth.(tstation[2],tstation[1],yrg',xrg) .- a0
drg = SOT.dist.(tstation[1], tstation[2], xrg, yrg')/1e6
aidx = abs.(drg.*sind.(arg)) .<= π*3.2*5/180
ntrg = length(trg)
τrg = zeros(ntrg)
for i = 1:ntrg
    τrg[i] = sum(aidx.*Δsrg[:,:,i])/sum(aidx)*3.2
end

K = -6
rc("font", size=8)
rc("axes", titlesize="medium")
fig,axs=plt.subplots(2,1,figsize=(190/25.4,190/25.4/2),sharex=true)
ax = axs[1]
p1, = ax.plot(te,τa/K, linewidth=1,color="tab:blue")
ax.plot(trg,τrg/K, linewidth=1, linestyle="--",color="tab:blue")
ax.fill_between(te, (τa-2ea)/K, (τa+2ea)/K, alpha=.2, zorder=3, color="tab:blue", linewidth=0)
p2, = ax.plot(te,τt/K, linewidth=1,color="tab:orange")
ax.fill_between(te, (τt-2et)/K, (τt+2et)/K, alpha=.2, zorder=3, color="tab:orange", linewidth=0)
p3, = ax.plot(te,τat/K, linewidth=1,color="tab:green")
ax.fill_between(te, (τat-2eat)/K, (τat+2eat)/K, alpha=.2, zorder=3, color="tab:green", linewidth=0)
ax.legend([p1, p2, p3], ["point mode 1 ","\$T\$ waves", "point with \$T\$ waves"], ncol=3, loc="upper center", frameon=false)
ax.set_ylabel("temperature anomaly (K)")
ax.set_title("\$19\\%\$(\$11\\%\$) full(trend) RMSE reduction for point with \$T\$ waves, \${10^\\circ}\$ subbasin")
ax = axs[2]
ax.plot(te,2ea/abs(K),color="tab:blue", linewidth=1)
ax.plot(te,2et/abs(K),color="tab:orange", linewidth=1)
ax.plot(te,2eat/abs(K),color="tab:green", linewidth=1)
ax.set_ylabel("error (K)")
ax.set_xlim([DateTime(2007,7,1),DateTime(2022,2,1)])
fig.align_ylabels()
fig.tight_layout()
fig.savefig("results/H11_argo_2d_tc.pdf")

####################
Random.seed!(3) # Setting the seed
d = Normal()
m = size(floats,1)
ns = σn*rand(d,m)

trd = Dates.value.(floats.t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
dx = SOT.dist.(floats.x', floats.y', floats.x, floats.y)/1e3

tm = trd[1]+(trd[m]-trd[1])/2
E = [I (trd.-tm)]
ω = 2π/SOT.meanyear
E = [E cos.(ω*trd) sin.(ω*trd) cos.(2ω*trd) sin.(2ω*trd)]

# solution covariance in time
C0 = σs^2*sqrt(2)*exp.(-abs.(trd.-trd')/λt-dx/λx/sqrt(2))
C = C0.*cos.(dx/λx/sqrt(2).-π/4) 

Cs = cholesky(sparse(C), perm=1:size(C,1))
Ls = sparse(Cs.L);
ss = Ls*rand(d, size(C,1))
y = E*[ss; zeros(5)]+ns

## Argo mean:
# 1112, 5481 profiles, mode 1
MLE: [57.31268805655394, 59.65000739540761, 0.4120606079727086, 0.023470327544687593]
lb: [52.41085872503784, 57.639992529373586, 0.39809427356708277, 0.019293569556395804]
ub: [62.672970681525264, 61.7301152573701, 0.42651692304291916, 0.028551288730929145]
# 1112, 5481 profiles, mode 2
MLE: [56.802200601029845, 59.35488252756207, 0.32463911758155617, 0.015638545849479005]
lb: [51.74782117001601, 57.21620252089478, 0.31312160196533406, 0.012393881121755686]
ub: [62.3502578498734, 61.57350408870858, 0.3365802806406162, 0.019732649836132432]

# small font
rc("font", size=8)
rc("axes", titlesize="medium")
N = 50

fig, ax = subplots(2, 2, figsize=(6.8, 6.8))
ax = reshape(ax,(4,))
xlabels = ["\$\\lambda_t\$ (days)","\$\\lambda_x\$ (km)","\$\\sigma_s\$ (ms/km)","\$\\sigma_\\eta\$ (ms)"]
db = [15,5,4e-5,1.5e-5]
k = 1
for i in [1,3,2,4]
    @printf("parameter %d...\n",k)
    xg = xk[:]
    gridi = LinRange(exp(xk[k])-db[k], exp(xk[k])+db[k], N)

    loglikei = similar(gridi)
    for (j,g) in enumerate(gridi)
        xg[k] = log(g)
        loglikei[j] = SOT.loglikelihoodargo(xg, y, σtrend, σannual, σsemiannual, floats.t, floats.x, floats.y;grad=false)
    end
    wi = SOT.exp_and_normalise(loglikei)
    if k>2
        gridi *= 1e3
    end
    ax[i].plot(gridi,wi)
    ax[i].set_xlabel(xlabels[k])
    ax[i].set_title("($(('a':'z')[k]))", loc="left")
    #ax[i].set_ylim([0,0.12])
    global k += 1
end
ax[1].set_ylabel("likelihood")
ax[2].set_ylabel("likelihood")
fig.tight_layout()
fig.savefig("results/marginalike_argo.pdf")
