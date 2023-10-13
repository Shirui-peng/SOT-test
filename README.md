# SOT-test
The code in this package processes seismic and acoustic data to infer temperature change in the ocean. The principal functions to be performed are:

Obtain an earthquake catalog from ISC.
Enhance the catalog through template matching at reference stations.
Find P-wave pairs at the reference stations.
Cross-correlate the waveforms of the identified P-wave pairs at the T-wave station.
Invert the measured travel time changes for travel time anomalies and apply a cycle skipping correction.
