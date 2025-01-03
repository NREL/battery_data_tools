![logo](assets/logo.png)

# battery_data_tools
This is a collection of battery data analysis tools to help users analyze battery test data, associated with the DOE EERE battery data hub, [batterydata.energy.gov](https://batterydata.energy.gov/).

## Handling battery data
Tools for handling battery data:
- `ampworks`
- `battery-data-toolkit`
- [batterydata.energy.gov formatting and plotting code](batterydata.energy.gov)

## DC pulse analysis
Methods for DC-pulse analysis
- Without a battery model: using simple analytical expressions using lmfit (but these neglect change of OCV due to charging/discharging).
- Make a simple model where OCV is assumed to just be linear, and is fit from the data, assuming we start at fully rested state
- Using a ECM+OCV model: DC pulses can be fit with an equivalent circuit model using the `thevenin` package and scipy. These will account for OCV change, and potentially hystersis, but require knowledge of the open-circuit-voltage curve of the battery.

## EIS analysis
EIS data can be plotted, fit, and analyzed using either a defined equivalent circuit model using the `impedance` package, or using the model agnostic distribution of relaxation times - distribution of phasances (DRT-DOP) method from the `hybrid-drt` package. See the [EIS Analysis Example](examples/EIS%20-%20ECM%20and%20DRT-DOP/Example%20EIS%20Analysis.ipynb).

## Differential voltage-capacity analysis
- Without a battery model: peak fitting using lmfit
- With a battery model: half-cell potential curve fitting using Ampworks. See the [DVQ Analysis Example](examples/Differential%20voltage-capacity/Example%20DVQ%20Analysis.ipynb)
- Simulation: [`alawa](https://www.hnei.hawaii.edu/alawa/) 

## Degradation trajectory prediction
Fitting of empirical models: `lmfit`, bayesian fitting

