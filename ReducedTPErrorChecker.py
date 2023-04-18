###############################################################################
# File: ReducedTPError.py
# Description: Tool used to check the expected time-integrated activity error
#              for a given sampling schedule based on simulated TIA data.
# Author: Avery B Peterson
# Institution: Wayne State University & University of Michigan
# Date: 17 APR 2023
# Python Version: 3.9.6
###############################################################################

# Imports
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import MultipleLocator, NullLocator


def add_gridlines(fig: plt.figure, alpha: float=1, which: str='both',
                  axis: str='both', behind: bool=False, fig_is_ax: bool=False):
    """
    Adds gridlines to a single axes or all axes in a matplotlib figure

    Parameters
    ----------
    fig : plt.figure
        Figure to add gridlines. Can be axes object if fig_is_ax==True
    alpha : float, optional
        Opacity of gridlines, by default 1
    which : str, optional
        Specifies which tick locations to draw gridlines. Valid options are
        'major', 'minor', or 'both', by default 'both'
    axis : str, optional
        Specifies which axis to draw gridlines. Valid options are 'x', 'y', or
        'both', by default 'both'
    behind : bool, optional
        If True, gridlines will be plotted behind axes elements, by default False
    fig_is_ax : bool, optional
        Specify True if fig is an axes object and not figure, by default False
    """
    def _add_gridlines_to_ax(ax, alpha, which, axis, behind):
        ax.set_axisbelow(behind)
        ax.minorticks_on()

        major_y_ticks = ax.yaxis.get_majorticklocs()
        major_y_spacing = major_y_ticks[1] - major_y_ticks[0]
        major_x_ticks = ax.xaxis.get_majorticklocs()
        major_x_spacing = major_x_ticks[1] - major_x_ticks[0]

        if (axis == 'both' or axis =='y') and (which == 'both' or which == 'minor'):
            ax.yaxis.set_minor_locator(MultipleLocator(major_y_spacing/2))
        else:
            ax.yaxis.set_minor_locator(NullLocator())

        if (axis == 'both' or axis == 'x') and (which == 'both' or which == 'minor'):
            ax.xaxis.set_minor_locator(MultipleLocator(major_x_spacing/2))
        else:
            ax.xaxis.set_minor_locator(NullLocator())

        ax.grid(True, which=which, axis=axis, alpha=alpha) 

    if fig_is_ax:
        _add_gridlines_to_ax(fig, alpha, which, axis, behind)
        return None
    for ax in fig.axes:
        _add_gridlines_to_ax(ax, alpha, which, axis, behind)
    return None


def find_ind_match(TPs: list):
    """
    Find sampling schedule in simulation data that most closely matches input
    sampling schedule and retrieve TIA error from appropriate CSV

    Parameters
    ----------
    TPs : list
        List of input imaging time points as floats

    Returns
    -------
    tuple
        Returns a tuple indicating the TIA error for that sampling schedule and
        the imaging time points in hours
    """
    # If single time point, read in, and return, both Hänscheid (H_diff.csv)
    # and Madsen (M_diff.csv) results.
    if len(TPs) == 1:
        df_H = pd.read_csv('H_diff.csv', index_col=0, header=[0])
        df_M = pd.read_csv('M_diff.csv', index_col=0, header=[0])
        # Find closest sampling time point
        df_inds_H = np.array([np.array(x, dtype=int) for x in df_H.columns])
        df_inds_M = np.array([np.array(x, dtype=int) for x in df_M.columns])
        match_H_ind = df_inds_H[np.abs(df_inds_H - TPs).argmin()]
        match_M_ind = df_inds_M[np.abs(df_inds_M - TPs).argmin()]

        # Set column labels to be int instead of str
        df_H.columns = df_inds_H.copy()
        df_M.columns = df_inds_M.copy()

        return (df_H.loc[:,match_H_ind], df_M.loc[:,match_M_ind]),\
               (match_H_ind, match_M_ind)

    # If given 2 or 3 time points, read in correct CSV
    if len(TPs) == 2:
        df = pd.read_csv('2TP_diff.csv', index_col=0, header=[0,1])
    if len(TPs) == 3:
        df = pd.read_csv('3TP_diff.csv', index_col=0, header=[0,1,2])
    # Find closest sampling schedule
    df_inds = np.array([np.array(x, dtype=int) for x in df.columns])
    match_ind = df_inds[np.abs(df_inds - TPs).sum(axis=1).argmin()]

    # Set column labels to be int instead of str
    df.columns = pd.MultiIndex.from_arrays(df_inds.T).copy()

    return df.loc[:,tuple(x for x in match_ind)], match_ind


def plot_stats(series: pd.Series, match_inds: tuple):
    """
    Generates plots for Kidney, Healthy Liver, Spleen, and Tumor displaying the
    TIA percent error distribution for the given sampling schedule. Also
    creates a dataframe with Mean, SD, 95% CI, Median, Min, and Max percent
    error

    Parameters
    ----------
    series : pd.Series
        Pandas Series (or tuple of Series) of the percent difference in TIA
        compared to ground truth that results from refitting at the indicated
        sampling schedule
    match_inds : tuple
        Tuple of floats that indicate the current sampling schedule
    """
    # Generates percent error plot and records stats for results dataframe
    def _plot_once(series, ax, title):
        plot_df = pd.concat((series.index.to_series(), series.copy()), axis=1, ignore_index=True)
        plot_df.columns = ['ROI', 'Diff']
        sns.boxplot(data=plot_df, x='ROI', y='Diff', ax=ax, whis=np.inf)
        plt.draw()
        results_df = pd.DataFrame(columns=['Structure','Mean (%)','SD (%)',
                                           '95% Confidence Interval (%)','Median (%)',
                                           'Minimum (%)','Maximum (%)'])
        for i, label in enumerate(ax.get_xticklabels()):
            x, y = label.get_position()
            roi = label.get_text()
            y = plot_df.loc[plot_df['ROI'] == roi].copy()
            mean = y['Diff'].mean()
            std = y['Diff'].std(ddof=1)

            ax.plot(x, mean, marker='X', color='k', zorder=3, markersize=8)
            ax.plot(x, mean + 1.96*std, marker=7, color=sns.color_palette()[i%4], markersize=8)
            ax.plot(x, mean - 1.96*std, marker=6, color=sns.color_palette()[i%4], markersize=8)
            ax.set_title(title)

            results_row = [roi, f'{mean:.1f}', f'{std:.1f}',
                           f'{mean-1.96*std:.1f}-{mean+1.96*std:.1f}',
                           f'{y["Diff"].median():.1f}', f'{y["Diff"].min():.1f}',
                           f'{y["Diff"].max():.1f}']

            results_df.loc[i,:] = results_row

        ax.set_xlabel(None)
        ax.set_ylabel('Percent error (%)')

        return results_df

    # If series is a tuple of Series, generate two plots (Hänscheid and Madsen)
    if isinstance(series, tuple):
        fig, axs = plt.subplots(1,2,figsize=(12,6))
        results_df = []
        for i,d in enumerate(series):
            plot_title = f'{"Hänscheid" if i==0 else "Madsen"} ({match_inds[i]}h)'
            results_df.append(_plot_once(d, axs[i], plot_title))
        # Combine two results dataframes into one
        empty_df = pd.DataFrame(columns=results_df[0].columns)
        rows_H = results_df[0].iterrows()
        rows_M = results_df[1].iterrows()
        for (i,row_H), (j,row_M) in zip(rows_H, rows_M):
            empty_df.loc[i,:] = [f'{x}/{y}' for x,y in zip(row_H,row_M)]
        results_df = empty_df.copy()
    # If series is just a Series, then generate plot and results dataframe once
    else:
        fig, ax = plt.subplots(1,1,figsize=(6,6))
        plot_title = f'{len(match_inds)} TP ('+', '.join([f'{x}h' for x in match_inds])+')'
        results_df = _plot_once(series, ax, plot_title)
    
    plt.draw()
    add_gridlines(fig, alpha=0.5, behind=True, axis='y')

    return results_df, fig


def create_table(df: pd.DataFrame, fontsize: int, title: str):
    """
    Creates a table as a matplotlib figure given an input dataframe

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame that contains data from table
    fontsize : int
        Specify the fontsize to use for the table
    title : str
        Provide text for the title of the table

    Returns
    -------
    plt.Figure
        Dataframe data as a table in a Matplotlib figure
    """
    fig, ax = plt.subplots(1,1,figsize=(14, 2))

    ax.axis('off')

    table = ax.table(cellText=df.values,
                     colLabels=df.columns,
                     loc='lower center')

    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)
    table.auto_set_column_width(range(len(df.columns)))
    table.scale(1,1.5)

    fig.text(0.5,0.8,title)

    return fig

def main():
    # Set figure params
    plt.rcdefaults()
    FONT_FAMILY = 'sans-serif'
    FONT = 'Arial'
    DEFAULT_SIZE = 14
    AXES_LABEL_SIZE = DEFAULT_SIZE
    TICK_LABEL_SIZE = DEFAULT_SIZE
    TITLE_SIZE = DEFAULT_SIZE
    MARKER_SIZE = 6
    MARKER_EDGE_WIDTH = 1
    GRID_WIDTH = 0.5
    LINE_WIDTH = 1.5
    BORDER_WIDTH = 2
    AXES_LABEL_PAD = 1
    MAJOR_TICK_SIZE = 3*BORDER_WIDTH

    rc('font', **{'family':FONT_FAMILY,
                  'sans-serif':[FONT],
                  'size':DEFAULT_SIZE})
    rc('figure', **{'titlesize':TITLE_SIZE})
    rc('axes', **{'titlesize':TITLE_SIZE,
                  'labelsize':AXES_LABEL_SIZE,
                  'labelpad':AXES_LABEL_PAD,
                  'linewidth':BORDER_WIDTH})
    rc('xtick', **{'labelsize':TICK_LABEL_SIZE,
                   'major.pad':AXES_LABEL_PAD,
                   'major.width':BORDER_WIDTH,
                   'major.size':MAJOR_TICK_SIZE,
                   'minor.width':BORDER_WIDTH/2,
                   'minor.size':MAJOR_TICK_SIZE/2})
    rc('ytick', **{'labelsize':TICK_LABEL_SIZE,
                   'major.pad':AXES_LABEL_PAD,
                   'major.width':BORDER_WIDTH,
                   'major.size':MAJOR_TICK_SIZE,
                   'minor.width':BORDER_WIDTH/2,
                   'minor.size':MAJOR_TICK_SIZE/2})
    rc('lines', **{'markersize':MARKER_SIZE,
                   'markeredgewidth':MARKER_EDGE_WIDTH,
                   'linewidth':LINE_WIDTH})
    rc('grid', **{'linewidth':GRID_WIDTH})

    # Loop until user manually stops
    while True:
        # Ask for and format sampling schedule
        TPs = input('\033[1mEnter up to three imaging times (in h) as a comma-separated list:\033[0m\n')
        # If input is 'exit' then exit program
        if TPs.lower() == 'exit':
            break
        # Remove spaces
        TPs.replace(' ','')
        # If text cannot be converted to float, ask user to try again
        try:
            TPs = [float(x) for x in TPs.split(',')]
        except ValueError:
            print('Input values could not be converted to float. '
                  'Please try again with different input or type \'exit\' to exit.\n')
            continue
        # Sort TPs
        TPs = sorted(TPs)

        # Find sampling schedule that is closest to input
        match_series, match_inds = find_ind_match(TPs)

        # Create stats figure and results dataframe for closest sampling schedule
        results_df, plot_fig = plot_stats(match_series, match_inds)

        # Notify user of closest match and generate title for table
        if len(match_inds) == 2 and match_inds[0] == match_inds[1]:
            match_string = f'{match_inds[0]}h'
            print(f'Closest match: {match_string}\n')
            table_title = f'Hänscheid/Madsen ({match_string})'
        else:
            match_string = ', '.join([f'{x}h' for x in match_inds])
            print(f'Closest match: {match_string}\n')
            table_title = f'{len(match_inds)} TP ({match_string})'
        
        # Create table figure
        table_fig = create_table(results_df, fontsize=DEFAULT_SIZE, title=table_title)

        # Display plots and table
        plt.show()

if __name__ == '__main__':
    # Run program and exit on KeyboardInterrupt
    try:
        main()
    except KeyboardInterrupt:
        pass
    print('Exiting...')