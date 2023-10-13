# SOT-test.jl
The code in this package processes seismic and acoustic data to infer temperature change in the Kuroshio Extension region. The principal functions to be performed are:

1. Obtain an earthquake catalog from ISC.
3. Find P-wave pairs at the reference stations.
4. Cross-correlate the waveforms of the identified P-wave pairs at the T-wave station.
5. Invert the measured travel time changes for travel time anomalies and apply a cycle skipping correction.
