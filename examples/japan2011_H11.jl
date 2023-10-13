include("../src/SOT.jl")
using .SOT, PyPlot, Printf, Dates, LinearAlgebra, Statistics, SparseArrays, DataFrames, CSV
using HDF5, Interpolations, Random, Distributions,NCDatasets

# identifier for experiment
eqname = "japan2011"
tsname = "H11"
evtpos = [38.10, 142.85]

# P-wave (reference) stations
pstations = ["IU.MAJO.00.BHZ","IU.MAJO.10.BHZ","PS.TSK..BHZ", "II.ERM.00.BHZ", "G.INU.00.BHZ"] 
pstnlats,pstnlons = [36.55,36.55,36.21,42.02,35.35],[138.2,138.2,140.11,143.16,137.03]
pstalocs = DataFrame(station=pstations,slat=pstnlats,slon=pstnlons)

# P-wave download source
psrc = "IRIS"

# intervals to which to cut P waveforms
pintervals = [[-3, 47], [-3, 47], [-3, 47], [-3, 47], [-3, 47]]

# frequency bands to which to filter P waveforms
pfreqbands = [[1.5, 2.5], [1.5, 2.5], [1, 3], [1.5, 2.5], [1.5, 2.5]]

# T-wave station
tstations = ["H11N3"]
tstalocs = [19.71786, 166.90986] 

# T-wave download source
tsrc = "IRIS"

# T-wave time window around predicted arrival time
tintervals = [[-20, 40]]

# T-wave filtering window width
tavgwidth = 0.5

# T-wave reference frequency at which to find first max CC
treffreq = 2.5

# frequencies used in inversion
tinvfreq = [2.5, 3.5]

# minimum CCs for T-wave pairs (at inversion frequencies)
tmincc = [0.6, 0.3]

# download P-wave data
SOT.downloadseisdata(eqname, tsname, pstations; src=psrc)

# cut and filter P waveforms
SOT.cutpwaves(eqname, tsname, pstations, pintervals, pfreqbands)

# find P-wave pairs
SOT.findpairs(eqname, tsname, pstations, pintervals, pfreqbands)
 
# measure T-wave lags Δτ
SOT.twavepickold(eqname, tsname, tstations, tintervals, tavgwidth, treffreq, pstations, pintervals, pfreqbands;soundspeed=1.47e3)

# collect usable pairs
tpairs, ppairs = SOT.collectpairs(eqname, tsname, tstations, tintervals, tavgwidth, treffreq,
                                  tinvfreq, tmincc, pstations, pintervals, pfreqbands;soundspeed=1.47e3)

# number of good T- and P-wave pairs
nt = size(tpairs, 1)
np = size(ppairs, 1)
    
# number of frequencies
l = length(tinvfreq)-1

# correlation time (days)
λt = 60

# correlation azimuth (degrees)
λθ = 2.0

# solution standard deviation for travel time anomalies (s)
στ = [0.27]

# location noise (s)
σx,σh = 0.019,0.028

# noise (s)
σn,σnp = 8.1e-3,2.2e-3

# origin time correction standard deviation (s)
σp = 0.974

# trend prior for coefficients of singular vectors (s/day)
σtrend = 0.01/SOT.meanyear

# annual cycle prior (s)
σannual = 0.1

# semi-annual cycle prior (s)
σsemiannual = 0.1

@printf("lt = %.0f days, n = %.2e s, nx = %.2e s\n\n",λt,σn,σx)

t, lon, lat, θ, E, R, N, P, D, invRyy = SOT.invertf1(tpairs, ppairs, tstalocs, pstalocs, evtpos, λt, λθ, στ, σx, σh, σn, σnp, σp; σtrend, σannual, σsemiannual)
m = length(t)

tpairs.Δτ,tpairs.Δτp,tpairs.cs,tpairs.x1,tpairs.x2,tpairs.Δτ2 = SOT.correctcycleskippingf1r2(eqname, tstations, tpairs, ppairs, E, R, N, P, m)

# collect delays into data vector
y = [reshape(vcat([(tpairs.Δτ[i])' for i = 1:nt]...), l*nt); ppairs.Δτ]
y2 = [tpairs.Δτ2; ppairs.Δτ]

tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24

# invert
a = P*E'*inv(Array(N))*y
a2 = P*E'*inv(Array(N))*y2
 
# extract trends
trends = a[2*m+1]
Ptrends = P[2*m+1,2*m+1]
@printf("\ntrend = %.2e K/yr, ptrend = %.2e\n",-trends*SOT.meanyear/6,2*sqrt(Ptrends)*SOT.meanyear/6)

# extract annual cycle
annual = a[2*m+2:2*m+3]
Pannual = diag(P[2*m+2:2*m+3,2*m+2:2*m+3])
@printf("\nannual = %s K, pannual = %s K\n",-annual/6,2*sqrt.(Pannual)/6)

# extract semiannual cycle
semiannual = a[2*m+4:2*m+5]
Psannual = diag(P[2*m+4:end,2*m+4:end])
@printf("\nsemiannual = %s K, psannual = %s K\n\n",-semiannual/6,2*sqrt.(Psannual)/6)

# reconstruct full travel time anomalies
τ = reshape(D*a, (m, 1)) 
eτ = reshape(sqrt.(diag(D*P*D')), (m, 1))

filename = @sprintf("results/pairs/%s_%s.h5",eqname,tsname)
# save to file
h5open(filename, "w") do file
  write(file, "t", Dates.value.(t .- DateTime(2000, 1, 1, 0, 0, 0)))
  write(file, "lon", lon)
  write(file, "lat", lat)
  write(file, "θ", θ)
  write(file, "τ", τ)
  write(file, "e", eτ)
  write(file, "y", y)
  write(file, "a", a)
end