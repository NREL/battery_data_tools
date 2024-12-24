# battery_data_tools
![logo](assets\logo.png)

This is a collection of battery data analysis tools to help users analyze battery test data, associated with the DOE EERE battery data hub, [batterydata.energy.gov](https://batterydata.energy.gov/).

## Handling battery data
`ampworks`, battery-data-toolkit, batterydata.energy.gov format definition

## DC pulse analysis
Without a battery model: using simple analytical expressions using lmfit (but these neglect change of OCV due to charging/discharging).
Make a simple model where OCV is assumed to just be linear, and is fit from the data, assuming we start at fully rested state?
Using a ECM+OCV model: DC pulses can be fit with an equivalent circuit model using the `thevenin` package and scipy. ALSO BLAST MATLAB CODE. These will account for OCV change, and potentially hystersis, but require knowledge of the open-circuit-voltage curve of the battery.

## EIS analysis
EIS data can be plotted, fit, and analyzed using either a defined equivalent circuit model using the `impedance` package, or using the model agnostic distribution of relaxation times - distribution of phasances (DRT-DOP) method from the `hybrid-drt` package.

## Differential voltage-capacity analysis
Without a battery model: peak fitting using lmfit
With a battery model: half-cell potential curve fitting using Ampworks
Simulation: [`alawa](https://www.hnei.hawaii.edu/alawa/) 

## Degradation trajectory prediction
Fitting of empirical models: `lmfit`, bayesian fitting (calendarAging repo)

