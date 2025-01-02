# Nina Prakash, NREL, 2024

import pandas as pd
import plotly.graph_objects as go
import plotly.colors as pcolors
from plotly.subplots import make_subplots
import numpy as np
from typing import Callable, Union
import random
import warnings


TEMPLATE='plotly_white'


parser = argparse.ArgumentParser()
parser.add_argument(
    "-p",
    "--data_processed",
    type=str,
    required=True,
    help="Path to the .xslx file with raw voltage-current data processed by process_iVt_data.m.",
)
parser.add_argument(
    "-s",
    "--data_summary",
    type=str,
    required=True,
    help="Path to the .xslx file with summarized battery performance metrics, calculated using process_IVt_data.m.",
)

#--------------------------------------------------------------------------------------------------------
def get_tick_label(tick_var):
    """
    Based on axis marker, function returns pre-defined 
    labels based on the passed in text
    
    Parameters
    -------------------------
    tick_var : str
         Thename of a particular column form a summary data or raw file
    
    Returns
    ------------------------
    tick_var_str : str
        Pre-defined, formatted string to be used for plotting.
    """
    
    if tick_var == 'Cycle_Index':
        tick_var_str = 'Cycle'
    elif tick_var == 'Charge_Throughput_Ah':
        tick_var_str = 'Charge throughput (Ah)'
    elif tick_var == 'Energy_Throughput_Wh':
        tick_var_str = 'Energy throughput (Wh)'
    elif tick_var == 'Absolute_Charge_Throughput_Ah':
        tick_var_str = 'Absolute Charge throughput (Ah)'
    elif tick_var == 'Absolute_Energy_Throughput_Wh':
        tick_var_str = 'Absolute Energy throughput (Wh)'
    elif tick_var == 'Equivalent_Full_Cycles':
        tick_var_str = 'Equiv. full cycles'
    elif tick_var == 'Time_s':
        tick_var_str = 'Time (s)'
    elif tick_var == 'tsecs_start':
        tick_var_str = 'Time (s)'
    elif tick_var == 'tsecs_end':
        tick_var_str = 'Time (s)'
    elif tick_var == 'Time_d':
        tick_var_str = 'Time (days)'
    elif tick_var ==  'Datenum_d':
        tick_var_str = 'Test datenum (days)'
    elif tick_var ==  'datenum_d':
        tick_var_str = 'Relative test date (days)'
    elif tick_var ==  'tsecs_cycle':
        tick_var_str = 'Cycle duration (s)'
    elif tick_var ==  'Q_chg':
        tick_var_str = 'Charge capacity (Ah)'
    elif tick_var ==  'q_chg':
        tick_var_str = 'Relative charge capacity'
    elif tick_var ==  'Q_dis':
        tick_var_str = 'Discharge capacity (Ah)'
    elif tick_var ==  'q_dis':
        tick_var_str = 'Relative discharge capacity'
    elif tick_var ==  'CE':
        tick_var_str = 'Coulombic efficieny'
    elif tick_var ==  'E_chg':
        tick_var_str = 'Charge energy (Wh)'
    elif tick_var ==  'e_chg':
        tick_var_str = 'Relative charge energy'
    elif tick_var ==  'E_dis':
        tick_var_str = 'Discharge energy (Wh)'
    elif tick_var ==  'e_dis':
        tick_var_str = 'Relative discharge energy'
    elif tick_var ==  'EE':
        tick_var_str = 'Energy efficiency'
    elif tick_var ==  'DeltaV':
        tick_var_str = 'Voltage hysteresis - dV (V)'
    elif tick_var ==  'deltaV':
        tick_var_str = 'Relative voltage hysteresis (V)'
    elif tick_var ==  'V_min':
        tick_var_str = 'Minimum voltage (V)'
    elif tick_var ==  'V_max':
        tick_var_str = 'Maximum voltage (V)'
    elif tick_var ==  'V_avg':
        tick_var_str = 'Average voltage (V)'
    elif tick_var ==  'I_min':
        tick_var_str = 'Minimum current (A)'
    elif tick_var ==  'I_max':
        tick_var_str = 'Maximum current (A)'
    elif tick_var ==  'I_avg':
        tick_var_str = 'Average current (A)'
    elif tick_var ==  'P_min':
        tick_var_str = 'Minimum power (W)'
    elif tick_var ==  'P_max':
        tick_var_str = 'Maximum power (W)'
    elif tick_var ==  'P_avg':
        tick_var_str = 'Average power (W)'
    elif tick_var ==  'T_min':
        tick_var_str = 'Minimum Temperature (^{\\circ}C)'
    elif tick_var ==  'T_max':
        tick_var_str = 'Maximum Temperature (^{\\circ}C)'
    elif tick_var ==  'T_avg':
        tick_var_str = 'Average Temperature (^{\\circ}C)'
    else:
        tick_var_str = tick_var

    return tick_var_str


#--------------------------------------------------------------------------------------------------------
def colormap(n, legend_style, plotly_colorscale = pcolors.sequential.Viridis):
    """
    Get an array of discrete colors of length n using Plotly builtin continuous colorscale 'Viridis'.
    If the legend_style is legend, randomize the color order so they're more distinguishable.
    
    Parameters
    -------------------------
    intermed : float
         Float in range [0, 1]
    
   Returns
   ------------------------
   colormap_array : tuple
        Color in RGB string tuple format
    """  
    intermed_values = np.linspace(0, 1, num=n)
    colors, _ = pcolors.convert_colors_to_same_type(plotly_colorscale)
    colorscale = pcolors.make_colorscale(colors)
    
    #Sub method to pull color
    def _get_intermediate_color(colorscale, intermed):

        if intermed <= 0 or len(colorscale) == 1:
            return colorscale[0][1]
        if intermed >= 1:
            return colorscale[-1][1]

        for cutoff, color in colorscale:
            if intermed > cutoff:
                low_cutoff, low_color = cutoff, color
            else:
                high_cutoff, high_color = cutoff, color
                break

        return pcolors.find_intermediate_color(
            lowcolor=low_color, highcolor=high_color,
            intermed=((intermed - low_cutoff) / (high_cutoff - low_cutoff)),
            colortype="rgb")
    
    colormap_array = [_get_intermediate_color(colorscale, i) for i in intermed_values]
    if legend_style == 'legend':
        random.seed(42)
        random.shuffle(colormap_array)

    return colormap_array


#--------------------------------------------------------------------------------------------------------
def plot_measurements(data_processed,
                             data_summary,
                             measurement_type,
                             fig,
                             direction: str = 'Both',
                             segment_label: str = 'All',
                             cycles: str = 'All',   
                             color_var: str = 'Cycle_Index',
                             linewidth: int = 2,
                             linestyle: str = 'solid',
                             linestyle_secondary: str = 'solid',
                             colormap: Callable[[int], str] = colormap,
                             plotly_colorscale = pcolors.sequential.Viridis,
                             axis_limits: list = [],
                             legend_style: str = 'legend',
                             legend_location: list = [1, 1]
                             ):
    """    
    Plots raw current-voltage-time battery data for various types ofmeasurements, 
    either plotting all recorded measurements or from up to 5 specific cycle 
    numbers. Can plot to an input axis handle, and optionally outputs the plotted axis.

    Parameters
    -------------------------
    data_processed: dataframe
          dataframe of raw current-voltage-time timeseries
          data. Must have been processed using the 'process_IVt_data'
          function before using plot_measurements.
    data_summary : dataframe
         dataframe of summarized battery performance
         metrics, calculated using process_IVt_data. Any SummaryData 
         field may be used to denote each measurement as set by the
         optional input 'ColorVariable'. Cycles to plot can be set by
         any Cycle_Label in the SummaryData.
    measurement_type : str 
         type of measurement to plot: 'Voltage-capacity', 'Voltage-time', 
        'Differential capacity', 'Differential voltage', 'EIS Nyquist', 
        'EIS Bode Real-Imaginary','EIS Bode Mag-Phase', 'Charge slippage'
    direction : str 
         'Charge', 'Discharge', or 'Both', used for plotting measurements 
         with clear charge or discharge directions (default 'Both').
    cycles : str 
         a Cycle_Label from SummaryData, or a vector of cycle numbers 
         (default 'All').
    color_var : str
         variable from summary data table to name legend entries or put as 
         colorbar tickets
    linewidth : int
         Width of line. Default = 1
    linestyle : str 
         Plotly dash line style, must be one of ("solid", "dot", "dash", 
         "longdash", "dashdot", or "longdashdot") or a dash length list in px 
         (eg. "5px,10px,2px,2px")
    linestyle_secondary : str 
         Plotly dash line style for secondary y axis, must be one of ("solid",
        "dot", "dash", "longdash", "dashdot", or "longdashdot")
    colormap :  str 
         Function handle that generates an arbitrarily long set of RGB colors.
    axis_limits : list of floats
         Axis limits as a list with values [x_min, x_max, y_min, y_max]   
    legend_style : str 
         Sets behavior for denoting each cycle as 'legend' (in a legend), 
         'colorbar' (using a colorbar), or 'off' (none) (default 'legend')
    legend_location : str 
         Sets location of legend as a list with values [x, y] where x sets the 
         x position with respect to xref from "0" (left) to "1" (right), and y sets 
         the y position with respect to yref from "0" (bottom) to "1" (top). 
         NOTE: using top left of the legend as the default anchor.
                                    
   Returns
   ------------------------
    fig : matplotlib/plotly figure
            the plotly go.Figure that was plotted on     
    """
    # Check for allowed input values
    _allowed_measurement_types = {'Voltage-capacity', 'Voltage-time',
                'Differential capacity', 'Differential voltage',
                'EIS Nyquist', 'EIS Bode Real-Imaginary',
                'EIS Bode Mag-Phase', 'Charge slippage'}
    if measurement_type not in _allowed_measurement_types:
        raise ValueError(f"`measurement_type` must be one of {_allowed_measurement_types}. Received {measurement_type}.")
        
    _allowed_directions = {'Discharge', 'Charge', 'Both'}
    if direction not in _allowed_directions:
        raise ValueError(f"`direction` must be one of {_allowed_directions}. Received {direction}.") 
    
    _allowed_color_vars = {'Cycle_Index','Absolute_Charge_Throughput_Ah','Absolute_Energy_Throughput_Wh',
        'Equivalent_Full_Cycles','Time_s','Time_d','Datenum_d','datenum_d','tsecs_start',
        'tsecs_end','tsecs_cycle','Q_chg','q_chg','Q_dis','q_dis','CE',
        'E_chg','e_chg','E_dis','e_dis','EE','DeltaV','deltaV','V_min','V_max','V_avg','I_min','I_max','I_avg',
        'P_min','P_max','P_avg','T_min','T_max','T_avg'}
    if color_var not in _allowed_color_vars:
        raise ValueError(f"`color_var` must be one of {_allowed_color_vars}. Received {color_var}.") 
    
    _allowed_linestyles = {"solid","dot", "dash", "longdash", "dashdot", "longdashdot"}
    if linestyle not in _allowed_linestyles:
        raise ValueError(f"`linestyle` must be one of {_allowed_linestyles}. Received {linestyle}.")
    if linestyle_secondary not in _allowed_linestyles:
        raise ValueError(f"`linestyle_secondary` must be one of {_allowed_linestyles}. Received {linestyle}.")
    
    if axis_limits:
        if len(axis_limits) != 4:
            raise ValueError(f"`axis_limits` must be a list of length 4 that represents [xmin, xmax, ymin, ymax]. Received list of length {len(axis_limits)}")
        
    _allowed_legend_styles = {"legend", "colorbar", "off"}
    if legend_style not in _allowed_legend_styles:
        raise ValueError(f"`legend_style` must be one of {_allowed_legend_styles}. Received {legend_style}.")
    
    if len(legend_location) != 2:
        raise ValueError(f"`legend_location` must be a list of length 2 that represents [x, y]. Received list of length {len(legend_location)}.")

    if segment_label != 'All':
        valid_segment_labels = data_processed['Segment_Label'].unique()
        if segment_label not in valid_segment_labels:
            raise ValueError(f"`segment_label` must be one of {valid_segment_labels}. Received {segment_label}.")

    if cycles == 'All':
        cycles = data_processed['Cycle_Index'].unique()
    else:
        if type(cycles) == str:
            valid_cycle_labels = data_summary['Cycle_Label'].unique()
            if cycles not in valid_cycle_labels:
                raise ValueError(f"Cycle must be one of {valid_cycle_labels}. Received {cycles}.")
            cycles = data_summary[data_summary['Cycle_Label'] == cycles]['Cycle_Index']
    
        else:
            # Check for integer type
            if any(not isinstance(cycle, int) for cycle in cycles):
                raise ValueError(f"List of cycle numbers must have integer type.")
            
            # Check that requested cycle indices are available
            for cycle in cycles:
                if cycle not in data_summary['Cycle_Index']:
                    raise ValueError(f"Cycle {cycle} is not available in the data_summary.")
    
    # TODO be able to pass in an existing "ax" / figure?
    # For now just creating a new one because the figure type
    # depends on the type of measurement...
    if measurement_type in {'EIS Bode Real-Imaginary', 'EIS Bode Mag-Phase'}:
        fig = make_subplots(specs=[[{"secondary_y": True}]])
    else:
        fig = go.Figure()
    fig.update_layout(template=TEMPLATE)
    
    # If n is small and this is an impedance bode plot, use the legend
    if len(cycles) < 10 or measurement_type in {'EIS Bode Real-Imaginary', 'EIS Bode Mag-Phase'}:
        legend_style = 'legend'
    else:
        legend_style = 'colorbar'
    
    # Line styles
    colors = colormap(len(cycles), legend_style, plotly_colorscale)
    tick_var_str = get_tick_label(color_var)
    
    # Plot lines
    trace_list = []

    color_var_values = []  # Save values for colorbar later
    for c, cycle in enumerate(cycles):

        # Split into charge and discharge data
        data_cycle = data_processed[data_processed['Cycle_Index'] == cycle]

        if segment_label != 'All':
            data_cycle = data_cycle[data_cycle['Segment_Label'] == segment_label]

        data_discharge = data_cycle[data_cycle['Current_A'] < 0]
        data_charge = data_cycle[data_cycle['Current_A'] > 0]
        
        if color_var != 'Cycle_Index':
            this_color_var = data_summary[data_summary['Cycle_Index'] == cycle][color_var].values[0]
        else:
            this_color_var = cycle
        color_var_values.append(this_color_var)

        if measurement_type == 'Voltage-capacity':
            if direction == 'Discharge':
                if not data_discharge.empty:
                    trace_list.append(go.Scatter(x=data_discharge['Discharge_Capacity_Ah'],
                                                    y=data_discharge['Voltage_V'],
                                                    line=dict(
                                                        width=linewidth,
                                                        color=colors[c],
                                                        dash=linestyle,
                                                    ),
                                                    name=f"{tick_var_str} {this_color_var}"))
            elif direction == 'Charge':
                if not data_charge.empty:
                    trace_list.append(go.Scatter(x=data_charge['Charge_Capacity_Ah'],
                                                    y=data_charge['Voltage_V'],
                                                    line=dict(
                                                        width=linewidth,
                                                        color=colors[c],
                                                        dash=linestyle,
                                                    ),
                                                    name=f"{tick_var_str} {this_color_var}"))
            elif direction == 'Both':
                if (not data_discharge.empty) and (not data_charge.empty):
                    Q_chg = data_charge['Charge_Capacity_Ah'].max()
                    trace_list.append(go.Scatter(x=Q_chg - data_discharge['Discharge_Capacity_Ah'],
                                                    y=data_discharge['Voltage_V'],
                                                    line=dict(
                                                        width=linewidth,
                                                        color=colors[c],
                                                        dash=linestyle,
                                                    ),
                                                    name=f"{tick_var_str} {this_color_var}"))
                    trace_list.append(go.Scatter(x=data_charge['Charge_Capacity_Ah'],
                                                    y=data_charge['Voltage_V'],
                                                    line=dict(
                                                        width=linewidth,
                                                        color=colors[c],
                                                        dash=linestyle_secondary,
                                                    ),
                                                    showlegend=False))

        elif measurement_type == 'Voltage-time':
            trace_list.append(go.Scatter(x=data_cycle['Time_s'] - data_cycle['Time_s'].values[0],
                                            y=data_cycle['Voltage_V'],
                                            line=dict(
                                                            width=linewidth,
                                                            color=colors[c],
                                                            dash=linestyle,
                                                        ),
                                        name=f"{tick_var_str} {this_color_var}"))
        elif measurement_type == 'Differential capacity':
            if direction == 'Discharge':
                if not data_discharge.empty:
                    trace_list.append(go.Scatter(x=data_discharge['Voltage_V'],
                                                    y=data_discharge['Differential_Capacity_Ah_V'],
                                                    line=dict(
                                                        width=linewidth,
                                                        color=colors[c],
                                                        dash=linestyle,
                                                    ),
                                                name=f"{tick_var_str} {this_color_var}"))
            elif direction == 'Charge':
                if not data_charge.empty:
                    trace_list.append(go.Scatter(x=data_charge['Voltage_V'],
                                                    y=data_charge['Differential_Capacity_Ah_V'],
                                                    line=dict(
                                                        width=linewidth,
                                                        color=colors[c],
                                                        dash=linestyle,
                                                    ),
                                                name=f"{tick_var_str} {this_color_var}"))
            elif direction == 'Both':
                if (not data_discharge.empty) and (not data_charge.empty):
                    trace_list.append(go.Scatter(x=data_discharge['Voltage_V'],
                                                    y=data_discharge['Differential_Capacity_Ah_V'],
                                                    line=dict(
                                                        width=linewidth,
                                                        color=colors[c],
                                                        dash=linestyle,
                                                    ),
                                                name=f"{tick_var_str} {this_color_var}"))
                    trace_list.append(go.Scatter(x=data_charge['Voltage_V'],
                                                    y=data_charge['Differential_Capacity_Ah_V'],
                                                    line=dict(
                                                        width=linewidth,
                                                        color=colors[c],
                                                        dash=linestyle,
                                                    ),
                                                showlegend=False))

        elif measurement_type == 'Differential voltage':
            if direction == 'Discharge':
                if not data_discharge.empty:
                    trace_list.append(go.Scatter(x=data_discharge['Discharge_Capacity_Ah'],
                                                    y=data_discharge['Differential_Voltage_V_Ah'],
                                                    line=dict(
                                                        width=linewidth,
                                                        color=colors[c],
                                                        dash=linestyle,
                                                    ),
                                                name=f"{tick_var_str} {this_color_var}"))

                
            elif direction == 'Charge':
                if not data_charge.empty:
                    trace_list.append(go.Scatter(x=data_charge['Charge_Capacity_Ah'],
                                                    y=data_charge['Differential_Voltage_V_Ah'],
                                                    line=dict(
                                                        width=linewidth,
                                                        color=colors[c],
                                                        dash=linestyle,
                                                    ),
                                                name=f"{tick_var_str} {this_color_var}"))

                
            elif direction == 'Both':
                if (not data_discharge.empty) and (not data_charge.empty):
                    trace_list.append(go.Scatter(x=data_discharge['Discharge_Capacity_Ah'],
                                                    y=data_discharge['Differential_Voltage_V_Ah'],
                                                    line=dict(
                                                        width=linewidth,
                                                        color=colors[c],
                                                        dash=linestyle,
                                                    ),
                                                name=f"{tick_var_str} {this_color_var}"))
                    trace_list.append(go.Scatter(x=data_charge['Charge_Capacity_Ah'],
                                                    y=data_charge['Differential_Voltage_V_Ah'],
                                                    line=dict(
                                                        width=linewidth,
                                                        color=colors[c],
                                                        dash=linestyle,
                                                    ),
                                                showlegend=False))

        elif measurement_type == 'EIS Nyquist':
            mask = data_cycle['Frequency_Hz'].notna()
            if mask.sum() == 0:
                warnings.warn(f"Cycle {cycle} has 'Frequency_Hz' values that are all NaN. Skipping plotting.")
                continue
            else:
                trace_list.append(go.Scatter(x=data_cycle['Z_Real_Ohm'],
                                            y=data_cycle['Z_Imag_Ohm'],
                                            line=dict(
                                                width=linewidth,
                                                color=colors[c],
                                                dash=linestyle,
                                            ),
                                        name=f"{tick_var_str} {this_color_var}",
                                        ))
                # Markers on frequency decades
                logfreq = round(np.log10(data_cycle['Frequency_Hz'][mask]))
                freq_decades = np.logspace(logfreq.min(), logfreq.max(), int(logfreq.max() - logfreq.min())+1)
                Z_Real_decades = np.interp(freq_decades, data_cycle['Frequency_Hz'][mask], data_cycle['Z_Real_Ohm'][mask])
                Z_Imag_decades = np.interp(freq_decades, data_cycle['Frequency_Hz'][mask], data_cycle['Z_Imag_Ohm'][mask])
                trace_list.append(go.Scatter(x=Z_Real_decades, y=Z_Imag_decades,
                                            mode='markers',
                                            marker=dict(
                                                color=colors[c],
                                                size=linewidth*3,
                                            ),
                                            showlegend=False))
                fig.update_yaxes(autorange="reversed")
            
        elif measurement_type == 'EIS Bode Real-Imaginary':
            secondary_y = False
            trace_list.append([go.Scatter(x=data_cycle['Frequency_Hz'],
                                            y=data_cycle['Z_Real_Ohm'],
                                            line=dict(
                                                width=linewidth,
                                                color=colors[c],
                                                dash=linestyle,
                                            ),
                                        name=f"{tick_var_str} {this_color_var} (Real)"),
                                        secondary_y])

            fig.update_xaxes(type="log")
            fig.update_yaxes(autorange="reversed")

            secondary_y = True
            trace_list.append([go.Scatter(x=data_cycle['Frequency_Hz'],
                                            y=data_cycle['Z_Imag_Ohm'],
                                            line=dict(
                                                width=linewidth,
                                                color=colors[c],
                                                dash=linestyle_secondary,
                                            ),
                                        name=f"{tick_var_str} {this_color_var} (Imag)"),
                                        secondary_y])
            
        elif measurement_type == 'EIS Bode Mag-Phase':
            trace_list.append(go.Scatter(x=data_cycle['Frequency_Hz'],
                                            y=data_cycle['Z_Mag_Ohm'],
                                            line=dict(
                                                width=linewidth,
                                                color=colors[c],
                                                dash=linestyle,
                                            ),
                                        name=f"{tick_var_str} {this_color_var} (|Z|)"))
            fig.update_xaxes(type="log")
            trace_list.append(go.Scatter(x=data_cycle['Frequency_Hz'],
                                            y=data_cycle['Z_Phase_Degree'],
                                            line=dict(
                                                width=linewidth,
                                                color=colors[c],
                                                dash=linestyle_secondary,
                                            ),
                                        name=f"{tick_var_str} {this_color_var} (\\Theta)"))
            
        elif measurement_type == 'Charge slippage':
            trace_list.append(go.Scatter(x=data_cycle['Charge_Throughput_Ah'],
                                            y=data_cycle['Voltage_V'],
                                            line=dict(
                                                width=linewidth,
                                                color=colors[c],
                                                dash=linestyle,
                                            ),
                                        name=f"{tick_var_str} {this_color_var}"))
    
        if measurement_type == 'Voltage-capacity':
            fig.update_layout(xaxis_title='Capacity (Ah)', yaxis_title='Voltage (V)')
        elif measurement_type == 'Voltage-time':
            fig.update_layout(xaxis_title='Time(s)', yaxis_title='Voltage (V)')
        elif measurement_type == 'Differential capacity':
            fig.update_layout(xaxis_title='Voltage (V)', yaxis_title='Differential capacity (A\\h\\cdotV^{-1})')
        elif measurement_type == 'Differential voltage':
            fig.update_layout(xaxis_title='Capacity (Ah)', yaxis_title='Differential voltage (V\\cdotAh^{-1})')
        elif measurement_type == 'EIS Nyquist':
            fig.update_layout(xaxis_title='Z_{Real} (\\Omega)', yaxis_title='Z_{Imaginary} (\\Omega)')
        elif measurement_type == 'EIS Bode Real-Imaginary':
            fig.update_xaxes(title='Frequency (Hz)')
            fig.update_yaxes(title='Z_{Real} (\\Omega)', secondary_y=False)
            fig.update_yaxes(title='Z_{Imaginary} (\\Omega)', secondary_y=True)
        elif measurement_type == 'EIS Bode Mag-Phase':
            fig.update_xaxes(title='Frequency (Hz)')
            fig.update_yaxes(title='|Z| (\\Omega)', secondary_y=False)
            fig.update_yaxes(title='\\Theta (\\circ)', secondary_y=True)
        elif measurement_type == 'Charge slippage':
            fig.update_layout(xaxis_title='Charge throughput (Ah)', yaxis_title='Voltage (V)')
                
    if axis_limits:
        xmin, xmax, ymin, ymax = axis_limits
        fig.update_layout(yaxis_range=[ymin, ymax], xaxis_range=[xmin, xmax])

    if legend_style == 'off':
        fig['layout']['showlegend'] = False

    if legend_style == 'colorbar':
        # Remove the default legend
        fig['layout']['showlegend'] = False
        # Add a manually-constructed colorbar
        cmin, cmax = min(color_var_values), max(color_var_values)
        cq1 = int((cmax-cmin)*.25)
        cmid = int((cmax-cmin)*.5)
        cq3 = int((cmax-cmin)*.75)
        x, y = legend_location
        colorbar_trace  = go.Scatter(x=[None],
                                y=[None],
                                mode='markers',
                                marker=dict(
                                    colorscale=plotly_colorscale, 
                                    showscale=True,
                                    cmin=cmin,
                                    cmax=cmax,
                                    colorbar=dict(thickness=20, tickvals=[cmin, cq1, cmid, cq3, cmax],
                                                  ticktext=[f"{tick_var_str} = {cmin}", 
                                                            f"{tick_var_str} = {cq1}",
                                                            f"{tick_var_str} = {cmid}",
                                                            f"{tick_var_str} = {cq3}",
                                                            f"{tick_var_str} = {cmax}"],
                                                  outlinewidth=0,
                                                  xanchor='left',
                                                  yanchor='top',
                                                  x=x,
                                                  y=y
                                                  )
                                ),
                                hoverinfo='none'
                                )
        fig.add_trace(colorbar_trace)

        if legend_style == 'legend':
            x, y = legend_location
            fig.update_layout(legend=dict(
                yanchor='top',
                xanchor='left',
                x=x,
                y=y,
            ))
    
        if measurement_type == 'EIS Bode Real-Imaginary':
            for trace in trace_list:
                fig.add_trace(trace[0], secondary_y=trace[1])
        else:
            for trace in trace_list:
                fig.add_trace(trace)
    
        #Provide more room for axis labels and set variable for autoscaling.
        fig.update_xaxes(title_standoff=10)
        fig.update_yaxes(title_standoff=20)
        fig.update_layout(autosize=True)
        return fig


#--------------------------------------------------------------------------------------------------------
def plot_metric(data_summary,
                x_var,
                y_var,
                legend_label = '',
                line_kwargs = {},
                cycles = 'All',
                axis_limits=[]
                ):
    """
    Plots the figures form the summary data files. Takes in the dataframe and
    associated inputs (e.g. axes) and returns the plotted figure. Note this 
    function works on one trace at a time. For multile trace plots see 
    plot_multiple_metrics().
    
    Useage: Q_figure = plotting.plot_metric(df, 'Cycle_Index', 'Q_dis')
    
    Parameters
    ---------------------------
    data_summary : df
          A single dataframe to plot
    x_var : str
         Variable name to plot on x-axis. Column from data_summary df.
    y_var : str
         Variable name to plot on y-axis. Column from data_summary df.
    legend_label : str
         Name of trace to display on legend
    line_kwargs : str
         Specifications to plotly line dict. Ex. linewidth, linestyle in
        ("solid", "dot", "dash", "longdash", "dashdot", or "longdashdot"). 
    cycles : str 
         'All', a Cycle_Label from data_summary, or a vector of cycle numbers (default 'All').  
    axis_limits :  list  of floats
         Axis limits as a list with values [x_min, x_max, y_min, y_max]      

    Returns
    ----------------------
       fig : matplotlib/plotly figure
            the plotly go.Figure that was plotted on     
    """

    if cycles != 'All':
        if isinstance(cycles, list):
            data_summary_cycles = data_summary[data_summary['Cycle_Index'].isin(cycles)]
        elif isinstance(cycles, str):
            if cycles not in data_summary['Cycle_Label'].unique():
                raise ValueError(f"Cycle must be one of {data_summary['Cycle_Label'].unique()}. Received {cycles}.")
            else:
                data_summary_cycles = data_summary[data_summary['Cycle_Label'] == cycles]
    else:
        data_summary_cycles = data_summary
    
    if x_var not in data_summary_cycles.columns:
        raise ValueError(f"x_var must be one of {data_summary_cycles.columns}. Received {x_var}.")
    if y_var not in data_summary_cycles.columns:
        raise ValueError(f"y_var must be one of {data_summary_cycles.columns}. Received {y_var}.")
        
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=data_summary_cycles[x_var],
                             y=data_summary_cycles[y_var],
                             name=legend_label,
                            line=dict(line_kwargs)))
    fig.update_layout(template=TEMPLATE,
                      title=f"{x_var} vs. {y_var}",
                      xaxis_title=get_tick_label(x_var),
                      yaxis_title=get_tick_label(y_var))
    fig.update_yaxes(title_standoff=20)
    if axis_limits:
        xmin, xmax, ymin, ymax = axis_limits
        fig.update_layout(yaxis_range=[ymin, ymax], xaxis_range=[xmin, xmax])
    return fig


# Reading excel is very slow- resaving as pickle for testing
# data_processed = pd.read_excel("data/processed data/P462 4 raw.xlsx")
# data_summary = pd.read_excel("data/processed data/P462 4 summary.xlsx")
# data_processed.to_pickle("data/processed data/P462 4 raw.pkl")
# data_summary.to_pickle("data/processed data/P462 4 summary.pkl")
#data_processed = pd.read_pickle("data/processed data/P462 2 raw.pkl")
#data_summary = pd.read_pickle("data/processed data/P462 2 summary.pkl")
#fig = plot_metric(data_summary, 'Time_s', 'Charge_Throughput_Ah', 'charge throughput',
                  #width=2, dash='longdashdot')
#fig = plot_measurements(data_processed, data_summary, 'Charge slippage',
                               #cycles=[1, 2, 3],
                               #linestyle='solid',
                               #legend_style='legend')
#fig.show()


#-------------------------------------------------------------------------------------------------------------
def plot_mutiple_metric(data_summary_list,
                x_var,
                y_var,
                legend_label = [],
                line_kwargs = {},
                cycles = 'All',
                axis_limits=[],
                ):
    """
    Method takes in a list of dataframes and plots each as a separate trace on
    the graph. Returns ocmpleted grapgh/plot object to the main function.
    
    Parameters
    ---------------------------
        data_summary_list : list
              A list of data frames to plot
        x_var : str
             Variable name to plot on x-axis. Column from data_summary df.
        y_var : str
             Variable name to plot on y-axis. Column from data_summary df.
        legend_label : list
             List of names of trace to display on legend
        line_kwargs : str
             Specifications to plotly line dict. Ex. linewidth, linestyle in
            ("solid", "dot", "dash", "longdash", "dashdot", or "longdashdot"). 
        cycles : str 
             'All', a Cycle_Label from data_summary, or a vector of cycle numbers (default 'All').  
        axis_limits :  list  of floats
             Axis limits as a list with values [x_min, x_max, y_min, y_max]      
        
    Returns
    ----------------------
       fig : matplotlib/plotly figure
            the plotly go.Figure that was plotted on     
    """
    fig = go.Figure()
    
    #Loop through all data frames. Add a new trace for each dataframe passed in
    for i, data_summary in enumerate(data_summary_list):
        #Parse dataframe based on cycles, or show all cycles
        if cycles != 'All':
            if isinstance(cycles, list):
                data_summary_cycles = data_summary[data_summary['Cycle_Index'].isin(cycles)]
            elif isinstance(cycles, str):
                if cycles not in data_summary['Cycle_Label'].unique():
                    raise ValueError(f"Cycle must be one of {data_summary['Cycle_Label'].unique()}. Received {cycles}.")
                else:
                    data_summary_cycles = data_summary[data_summary['Cycle_Label'] == cycles]
        else:
            data_summary_cycles = data_summary 
        
        #Validate that the passed in X and Y selected columns are part of the dataframe.
        if x_var not in data_summary_cycles.columns:
            raise ValueError(f"x_var must be one of {data_summary_cycles.columns}. Received {x_var}.")
        if y_var not in data_summary_cycles.columns:
            raise ValueError(f"y_var must be one of {data_summary_cycles.columns}. Received {y_var}.")
        
        #Build trace to add to figure.    
        fig.add_trace(go.Scatter(x=data_summary_cycles[x_var],
                                 y=data_summary_cycles[y_var],
                                 name=legend_label[i],
                                 line=dict(line_kwargs)
                                ))
    #Update layout with the Title, and X and Y axis lable
    fig.update_layout(template=TEMPLATE,
                      title=f"{x_var} vs. {y_var}",
                       xaxis_title=get_tick_label(x_var),
                       yaxis_title=get_tick_label(y_var),                       
                      )
    #Stand off the X and Y labels from the axis values a few pixels
    fig.update_xaxes(title_standoff=10)
    fig.update_yaxes(title_standoff=20)
    
    #Set plot limits if axis_limit values are passed in.
    if axis_limits:
        xmin, xmax, ymin, ymax = axis_limits
        fig.update_layout(yaxis_range=[ymin, ymax], xaxis_range=[xmin, xmax])
        
    #fig.show()
    #return
    return fig


if __name__ == "__main__":

    args = parser.parse_args()
    data_processed = pd.read_pickle(args.data_processed)
    data_summary = pd.read_pickle(args.data_summary)

    fig = plot_metric(data_summary, 'Time_s', 'Charge_Throughput_Ah', 'charge throughput',
                    width=2, dash='longdashdot',)
    fig.show()


    fig = plot_measurements(data_processed, data_summary, 'Charge slippage',
                                cycles='All',
                                segment_label='Pulse 50SOC',
                                linestyle='solid',
                                legend_style='legend')
    fig.show()

    print()
