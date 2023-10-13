# TODO:
# - Allow for clock error corrections (before cutting).

# Earth's radius
const earthradius = 6371e3

# sound speed in water (for arrival time estimation)
# const soundspeed = 1.46e3

"""
    twavepickold(eqname, tstations, tintervals, tavgwidth, treffreq, pstations,
              pintervals, pfreqbands; saveplot=false)
Picks the lags at which the cross-correlation function between pairs is maximized, saving
the results to file.
# Arguments
- `eqname::String`: earthquake name to identify experiment.
- `tstations::Array{String,1}`: *T*-wave station designations.
- `tintervals::Array{Array{Float64,1}}`: beginning and end of waveform relative to predicted
  arrival time.
- `tavgwidth::Float64`: width of Gaussian frequency window.
- `treffreq::Float64`: reference frequency at which to first find max cross-correlation.
- `pstations`, `pintervals`, `pfreqbands`: *P*-wave stations and parameters to specify pair
  catalogs.
- `saveplot::Boolean`: whether to save a plot.
# Example
```
julia> twavepick("nias", ["H08S2"], [[-10, 70]], [0.5], [2.0], ["PS.PSI..BHZ", "MY.KUM..BHZ", "II.WRAB.00.BHZ"], [[-3, 47], [-3, 47], [-3, 47]], [[1, 3], [1, 3], [1.5, 2.5]])
[...]
```
"""
function twavepickold(eqname, tsname, tstations, tintervals, tavgwidth, treffreq, pstations,
                      pintervals, pfreqbands; soundspeed=1.5e3, saveplot=false)

  # TODO: make sure start time correction is right
  # TODO: deal with data gaps
  # TODO: check that max tracing is working as intended

  # frequencies at which to measure delays
  tfreq = 0.1:0.1:10

  # load and combine pairss of P-wave pairs
  pairs = DataFrame[]
  for i = 1 : size(pstations, 1)
    filename = paircatfileold(eqname, tsname, pstations[i], pintervals[i], pfreqbands[i])
    push!(pairs, DataFrame(CSV.File(filename, select=1:10, comment="#")))
  end
  pairs = sort(unique(vcat(pairs...)))
  
  @printf("\npair number: %d\n\n", length(pairs.event1))

  # loop over T-wave stations
  for i = 1 : length(tstations)
    mkpath(tdelaydir(eqname, tstations[i], tintervals[i], tavgwidth, treffreq, soundspeed))
    
    # iterate over pairs
    for j = 1 : size(pairs, 1)

      @printf("%d/%d: %s %s\n", j, size(pairs, 1), pairs[j,:event1], pairs[j,:event2])

      # file to which measurements are saved
      tdelayfile = @sprintf("%s/%s_%s.h5",
                            tdelaydir(eqname, tstations[i], tintervals[i], tavgwidth,
                                      treffreq, soundspeed),
                            fmttime(pairs[j,:event1]), fmttime(pairs[j,:event2]))

      # file to which plot is saved
      if saveplot
        mkpath(tplotdir(eqname, tstations[i], tintervals[i], tavgwidth, treffreq, soundspeed))
        tplotfile = @sprintf("%s/%s_%s.pdf",
                             tplotdir(eqname, tstations[i], tintervals[i], tavgwidth,
                                      treffreq, soundspeed),
                             fmttime(pairs[j,:event1]), fmttime(pairs[j,:event2]))
      end

      # files from which to read T-wave data
      tdatafile1 = tdatafile(eqname, tstations[i], pairs[j,:event1])
      tdatafile2 = tdatafile(eqname, tstations[i], pairs[j,:event2])
      if ishydr(tstations[i]) && pairs[j,:event2]>=DateTime(2018, 1, 1)
        tdatafile2 = tdatafile(eqname, "IM.H11N3..EDH", pairs[j,:event2])
        if pairs[j,:event1]>=DateTime(2018, 1, 1)
          tdatafile1 = tdatafile(eqname, "IM.H11N3..EDH", pairs[j,:event1])
        end
      end

      # check whether all data is present
      if !isfile(tdelayfile) && isfile(tdatafile1) && isfile(tdatafile2) && filesize(tdatafile1) > 0 && filesize(tdatafile2) > 0

        # station location
        if ishydr(tstations[i])
          stalat, stalon = DataFrame(CSV.File(tstationlocfile(tstations[i])))[1,:]
        else
          stalat, stalon = h5read(tdatafile1, "latitude"), h5read(tdatafile2, "longitude")
        end

        # estimate travel time from geodetic distance
        evtlat = .5(pairs[j,:latitude1] + pairs[j,:latitude2])
        evtlon = .5(pairs[j,:longitude1] + pairs[j,:longitude2])
        range = dist(evtlon, evtlat, stalon, stalat)
        traveltime = range/soundspeed

        if ishydr(tstations[i])
          if pairs[j,:event1]>=DateTime(2018, 1, 1)
            # open files
            fid1 = h5open(tdatafile1, "r")
            fid2 = h5open(tdatafile2, "r")
          
            try
   
              # sampling interval
              Δ1 = read(fid1, "fs")^-1
              Δ2 = read(fid2, "fs")^-1
    
              # extract signals
              s1 = read(fid1, "trace")
              s2 = read(fid2, "trace")
    
              # first times in files
              starttime1 = DateTime(1970, 1, 1) + Microsecond(read(fid1, "starttime"))
              starttime2 = DateTime(1970, 1, 1) + Microsecond(read(fid2, "starttime"))

              # close files
              close(fid1)
              close(fid2)
            catch y
              # close files
              close(fid1)
              close(fid2)
              continue
            end
          else
            # load both waveforms
            trace1 = read_sac(tdatafile1)
          
            # sampling interval
            Δ1 = trace1.delta
          
            # first times in files
            starttime1 = trace1.evt.time + Microsecond(Int(round(trace1.b*1e6)))
          
            # extract signals
            s1 = Float64.(trace1.t)
            
            if pairs[j,:event2]>=DateTime(2018, 1, 1)
              fid2 = h5open(tdatafile2, "r")
          
              try
                Δ2 = read(fid2, "fs")^-1
                s2 = read(fid2, "trace")
                starttime2 = DateTime(1970, 1, 1) + Microsecond(read(fid2, "starttime"))
                close(fid2)
              catch y
                close(fid2)
                continue
              end
            else
              trace2 = read_sac(tdatafile2)
              Δ2 = trace2.delta
              starttime2 = trace2.evt.time + Microsecond(Int(round(trace2.b*1e6)))
              s2 = Float64.(trace2.t)
            end
          end
        else
          if j==1
            @printf("T wave in land station...\n")
          end
          # open files
          fid1 = h5open(tdatafile1, "r")
          fid2 = h5open(tdatafile2, "r")
          
          try

              # sampling interval
              Δ1 = read(fid1, "fs")^-1
              Δ2 = read(fid2, "fs")^-1
    
              # extract signals
              s1 = read(fid1, "trace")
              s2 = read(fid2, "trace")
    
              # first times in files
              starttime1 = DateTime(1970, 1, 1) + Microsecond(read(fid1, "starttime"))
              starttime2 = DateTime(1970, 1, 1) + Microsecond(read(fid2, "starttime"))
              
              if starttime1 < DateTime(2014, 7, 13) && starttime2 > DateTime(2014, 7, 13) && occursin("WAKE.00", tstations[i])#abs(Δ1-2*Δ2)<1e-3
                @printf("Downsampling event2...\n")
                Δ2 = read(fid1, "fs")^-1
                s2 = s2[1:2:end]
              end
              
              # close files
              close(fid1)
              close(fid2)
          catch y
              println(y.msg)
              # close files
              close(fid1)
              close(fid2)
              #rm(tdatafile1)
              #rm(tdatafile2)
              continue
          end

        end

        # average sampling interval
        Δ = (Δ1 + Δ2)/2

        # lengths of traces
        N1 = length(s1)
        N2 = length(s2)

        # times of samples
        time1 = starttime1 .+ Microsecond.(Int.(round.((0:N1-1)*Δ1*1e6)))
        time2 = starttime2 .+ Microsecond.(Int.(round.((0:N2-1)*Δ2*1e6)))

        # convert time of samples to seconds since event time
        time1 = Dates.value.(time1 .- pairs[j,:event1])/1e3
        time2 = Dates.value.(time2 .- pairs[j,:event2])/1e3

        # find first index of window
        i1 = findfirst(traveltime + tintervals[i][1] .< time1)
        i2 = findfirst(traveltime + tintervals[i][1] .< time2)

        # window length
        n1 = Int(round((tintervals[i][2] - tintervals[i][1])/Δ1))
        n2 = Int(round((tintervals[i][2] - tintervals[i][1])/Δ2))

        # check whether time series are long enough and window has same length
        if !isnothing(i1) && !isnothing(i2) && i1 + n1 - 1 ≤ N1 && i2 + n2 - 1 ≤ N2 && n1 == n2

          # sample length
          n = n1

          # make sure sample length is even
          (mod(n, 2) == 1) && (n -= 1)

          # cut signals
          s1 = s1[i1:i1+n-1]
          s2 = s2[i2:i2+n-1]

          # remove means
          s1 .-= mean(s1)
          s2 .-= mean(s2)

          # Hann window
          hann = .5*(1 .- cos.(2π*(1:n)/n))

          # Fourier transform
          st1 = rfft(hann.*s1)
          st2 = rfft(hann.*s2)

          # sample frequencies
          ω = (0:n÷2)/(n*Δ)

          # Gaussian for freq. averages
          G = exp.(-(ω.-tfreq').^2/2tavgwidth^2)

          # filter
          st1f = st1.*G
          st2f = st2.*G

          # normalization
          norm1 = (abs.(st1f[1,:]).^2 + 2sum(abs.(st1f[2:n÷2,:]).^2, dims=1)[1,:] + abs.(st1f[n÷2+1,:]).^2)/n
          norm2 = (abs.(st2f[1,:]).^2 + 2sum(abs.(st2f[2:n÷2,:]).^2, dims=1)[1,:] + abs.(st2f[n÷2+1,:]).^2)/n

          # cross-correlation function
          cc = circshift(irfft(conj(st1f).*st2f, n, 1)./sqrt.(norm1.*norm2)', (n÷2, 0))

          # offsets
          lags = collect(-n÷2:n÷2-1)*Δ

          # correct for difference in initial times
          lags .+= time2[i2] - time1[i1]

          # pick max CC and adjacent max
          Δτl, Δτl2, Δτc, Δτr, Δτr2, ccc, ccr, ccr2, ccl, ccl2 = findmaxcc(cc, tfreq, lags, treffreq, 1/Δ)
          i0 = argmin(abs.(tfreq .- treffreq))

          @printf("CC = %4.2f\n", ccc[i0])

#          # experimental: maximize over phase and group delays
#          δ = -1:.001:1;
#          G = exp.(-(ω.-2).^2/2tavgwidth^2);
#          s = exp.(2π*im*(ω.-2)/2 .* δ'.*ω);
#          ccg = circshift(irfft(conj.(st1).*st2.*G.^2 .* s, n, 1), (n÷2, 0))
#          Δτ = lags[argmax(ccg)[1]]
#          Δτg = δ[argmax(ccg)[2]] + Δτ
#          println("old: Δτ = $(Δτc[tfreq.==2][1])")
#          println("new: Δτ = $Δτ, Δτg = $Δτg")

          # save figure if CC ≥ 0.6
          if saveplot && ccc[i0] ≥ 0.6
            fig = figure()
            imshow(cc', cmap="RdBu_r", vmin=-1, vmax=1, origin="lower", aspect="auto",
                   extent=[lags[1] - Δ/2, lags[end] + Δ/2, tfreq[1] - .5(tfreq[2]-tfreq[1]),
                           tfreq[end] + .5(tfreq[end]-tfreq[end-1])])
            title(@sprintf("%.2f: %.2f s, %.2f: %.2f s, %.2f: %.2f s", ccl[i0], Δτl[i0], ccc[i0], Δτc[i0], ccr[i0], Δτr[i0]))
            plot(Δτl, tfreq, color="black")
            plot(Δτc, tfreq, color="black")
            plot(Δτr, tfreq, color="black")
            xlim(Δτc[i0]-1, Δτc[i0]+1)
            plt.colorbar()
            xlabel("lag (s)")
            ylabel("frequency (Hz)")
            savefig(tplotfile, dpi=200)
            close(fig)
            
            fig, ax = subplots(1, 1)
            ax.plot(lags,cc[:,i0])
            ax.set_xlabel("lag (s)")
            ax.set_ylabel("2.5 Hz CC")
            ax.set_xlim(Δτc[i0]-1, Δτc[i0]+1)
            fig.savefig(string(tplotfile[1:end-4],"2.5HzCC.pdf"), dpi=200)
            close(fig)
            
            designmethod = Butterworth(4)
            responsetype = Bandpass(treffreq-tavgwidth, treffreq+tavgwidth; fs=Int(round(1/Δ)))
            s1plt = filtfilt(digitalfilter(responsetype, designmethod), s1)
            s2plt = filtfilt(digitalfilter(responsetype, designmethod), s2)
            tplt = (0:n-1)*Δ
            fig, ax = subplots(2, 1)
            ax[1].plot(tplt, s1plt, color="black", alpha=0.75, label=string("event1,M",pairs[j,:magnitude1]))
            ax[1].plot(tplt, s2plt, color="red", alpha=0.75, label=string("event2,M",pairs[j,:magnitude2]))
            ax[1].legend(frameon=false, loc=4)
            ax[2].plot(tplt, s1plt/max(s1plt...), color="black", alpha=0.75)
            ax[2].plot(tplt, s2plt/max(s2plt...), color="red", alpha=0.5)
            ax[2].plot(tplt.-Δτc[i0], s2plt/max(s2plt...), color="green", alpha=0.75, label="event2 shifted")
            ax[2].legend(frameon=false, loc=4)
            ax[2].set_xlabel("time (s)")
            ax[1].set_ylabel("absolute amplitude")
            ax[2].set_ylabel("normalized amplitude")
            fig.savefig(string(tplotfile[1:end-4],"waveform.pdf"), dpi=200)
            close(fig)
          end

          # save to file
          h5open(tdelayfile, "w") do fid
            write(fid, "freq", collect(tfreq))
            write(fid, "Δτl", Δτl)
            write(fid, "Δτl2", Δτl2)
            write(fid, "Δτc", Δτc)
            write(fid, "Δτr", Δτr)
            write(fid, "Δτr2", Δτr2)
            write(fid, "ccc", ccc)
            write(fid, "ccr", ccr)
            write(fid, "ccr2", ccr2)
            write(fid, "ccl", ccl)
            write(fid, "ccl2", ccl)
#            # also save full info for plot
#            write(fid, "Δ", Δ)
#            write(fid, "cc", cc)
          end

        else

          # create empty file to prevent repeated measuring attempt
          touch(tdelayfile)

        end

      end

    end

  end

end

"Calculate great-circle distance"
function dist(lon1, lat1, lon2, lat2)
  cosΔσ = sind(lat1)*sind(lat2) + cosd(lat1)*cosd(lat2)*cosd(lon1-lon2)
  return earthradius*acos(sign(cosΔσ)*min(1,abs(cosΔσ)))
end

"Find maxima in CC starting at a reference frequency"
function findmaxcc(cc, freq, lags, reffreq, fs)

  # sample length
  n = length(lags)

  # number of central frequencies
  m = length(freq)

  # initialize lags, CCs
  Δτl = Array{Float64,1}(undef, m)
  Δτl2 = Array{Float64,1}(undef, m)
  Δτc = Array{Float64,1}(undef, m)
  Δτr = Array{Float64,1}(undef, m)
  Δτr2 = Array{Float64,1}(undef, m)
  ccc = Array{Float64,1}(undef, m)
  ccr = Array{Float64,1}(undef, m)
  ccr2 = Array{Float64,1}(undef, m)
  ccl = Array{Float64,1}(undef, m)
  ccl2 = Array{Float64,1}(undef, m)

  # frequency closest to the reference frequency
  i0 = argmin(abs.(freq.-reffreq))

  # index of max CC at reference frequency
  imaxc = argmax(cc[:,i0])

  # find max by interpolation
  Δτc[i0] = lags[imaxc] + argmaxquad(cc[mod1.(imaxc-1:imaxc+1, n),i0], 1/fs)
  ccc[i0] = maxquad(cc[mod1.(imaxc-1:imaxc+1, n),i0])

  # sweep frequencies
  sweep!(Δτc, ccc, cc, freq, lags, i0, imaxc, fs)

  # index of right secondary max CC at reference frequency
  idx = mod1.(imaxc .+ (Int(round(fs/2reffreq)) : Int(round(3fs/2reffreq))), n)
  imaxr = idx[argmax(cc[idx,i0])]

  # find max by interpolation
  Δτr[i0] = lags[imaxr] + argmaxquad(cc[mod1.(imaxr-1:imaxr+1, n),i0], 1/fs)
  ccr[i0] = maxquad(cc[mod1.(imaxr-1:imaxr+1, n),i0])

  # sweep frequencies
  sweep!(Δτr, ccr, cc, freq, lags, i0, imaxr, fs)
  
  # index of right secondary max CC at reference frequency
  idx = mod1.(imaxc .+ (Int(round(3fs/2reffreq)) : Int(round(5fs/2reffreq))), n)
  imaxr2 = idx[argmax(cc[idx,i0])]

  # find max by interpolation
  Δτr2[i0] = lags[imaxr2] + argmaxquad(cc[mod1.(imaxr2-1:imaxr2+1, n),i0], 1/fs)
  ccr2[i0] = maxquad(cc[mod1.(imaxr2-1:imaxr2+1, n),i0])

  # sweep frequencies
  sweep!(Δτr2, ccr2, cc, freq, lags, i0, imaxr2, fs)

  # index of left secondary max CC at reference frequency
  idx = mod1.(imaxc .+ (-Int(round(3fs/2reffreq)) : -Int(round(fs/2reffreq))), n)
  imaxl = idx[argmax(cc[idx,i0])]

  # find max by interpolation
  Δτl[i0] = lags[imaxl] + argmaxquad(cc[mod1.(imaxl-1:imaxl+1, n),i0], 1/fs)
  ccl[i0] = maxquad(cc[mod1.(imaxl-1:imaxl+1, n),i0])

  # sweep frequencies
  sweep!(Δτl, ccl, cc, freq, lags, i0, imaxl, fs)
  
  # index of left secondary max CC at reference frequency
  idx = mod1.(imaxc .+ (-Int(round(5fs/2reffreq)) : -Int(round(3fs/2reffreq))), n)
  imaxl2 = idx[argmax(cc[idx,i0])]

  # find max by interpolation
  Δτl2[i0] = lags[imaxl2] + argmaxquad(cc[mod1.(imaxl2-1:imaxl2+1, n),i0], 1/fs)
  ccl2[i0] = maxquad(cc[mod1.(imaxl2-1:imaxl2+1, n),i0])

  # sweep frequencies
  sweep!(Δτl2, ccl2, cc, freq, lags, i0, imaxl2, fs)

  return Δτl, Δτl2, Δτc, Δτr, Δτr2, ccc, ccr, ccr2, ccl, ccl2

end

"Get maximum of CC function by fitting a quadratic to three points"
maxquad(c) = c[2] + (c[1]-c[3])^2/(16c[2]-8(c[1]+c[3]))

"Get index of maximum of CC function by fitting a quadratic to three points"
argmaxquad(c, h) = h*(c[1]-c[3])/2(c[1]-2c[2]+c[3])

"Sweep across frequencies to trace max CC"
function sweep!(maxlags, ccs, cc, freq, lags, i0, imax, fs)

  # sample length
  n = length(lags)

  # number of central frequencies
  m = length(freq)

  # save reference index of maximum
  imaxref = imax

  # sweep up
  for i = i0+1:m

    # search region: ± a quarter period
    idx = mod1.(imax .+ (-Int(round(fs/4freq[i])) : Int(round(fs/4freq[i]))), n)

    # maximum on grid
    imax = idx[argmax(cc[idx,i])]

    # interpolated max
    maxlags[i] = lags[imax] + argmaxquad(cc[mod1.(imax-1:imax+1, n),i], 1/fs)
    ccs[i] = maxquad(cc[mod1.(imax-1:imax+1, n),i])

  end

  # restore reference index
  imax = imaxref

  # sweep down
  for i = i0-1:-1:1

    # search region: ± a quarter period
    idx = mod1.(imax .+ (-Int(round(fs/4freq[i])) : Int(round(fs/4freq[i]))), n)

    # maximum on grid
    imax = idx[argmax(cc[idx,i])]

    # interpolated max
    maxlags[i] = lags[imax] + argmaxquad(cc[mod1.(imax-1:imax+1, n),i], 1/fs)
    ccs[i] = maxquad(cc[mod1.(imax-1:imax+1, n),i])

  end

end

"Check if station is a CTBTO hydrophone station"
ishydr(station) = occursin(r"H[0,1][1-9][E,W,N,S][1-3]", station)

"Check if station is a 3-in-1 CTBTO hydrophone station"
arehydr(station) = occursin(r"H[0,1][1-9][E,W,N,S]", station)

function tdatafile(eqname, station, date)
  if ishydr(station) && date<DateTime(2018,1,1)
    return @sprintf("data/hydrdata/%s/%d/%d_%d.sac", station, year(date), year(date),
                    dayofyear(date))
  else
    fmttime = Dates.format(date, "yyyy-mm-ddTHH:MM:SS.ss")
    return @sprintf("%s/%s.h5", seisdatadir(eqname, station), fmttime)
  end
end

tstationlocfile(station) = @sprintf("data/hydrdata/%s/loc.csv", station)

tdelaydir(eqname, station, interval, avgwidth, reffreq, soundspeed) = @sprintf("data/tdelays/%s_%s_%+02d_%+02d_%3.1f_%3.1f_%.2f", eqname, station, interval[1], interval[2], avgwidth, reffreq, 1e-3soundspeed)

tplotdir(eqname, station, interval, avgwidth, reffreq, soundspeed) = @sprintf("data/tplots/%s_%s_%+02d_%+02d_%3.1f_%3.1f_%.2f", eqname, station, interval[1], interval[2], avgwidth, reffreq, 1e-3soundspeed)
