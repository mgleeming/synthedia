import os

from matplotlib.figure import Figure
from matplotlib.patches import Rectangle, Patch
import numpy as np
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
from .mzml import MZMLReader

def plot_tic(options):

    fig = Figure()
    ax = fig.subplots(2,1, sharex = True)

    for groupi in range(options.n_groups):
        for samplei in range(options.samples_per_group):

            mzml_file = os.path.join(
                options.out_dir, '%s_group_%s_sample_%s.mzML' %(
                    options.output_label, groupi, samplei
                )
            )

            ms1_rts, ms1_ints = [],[]
            ms2_rts, ms2_ints = [],[]

            for (rt, lvl, mzs, ints) in MZMLReader(mzml_file):

                if lvl == 1:
                    ms1_rts.append(rt)
                    ms1_ints.append(sum(ints))
                if lvl == 2:
                    ms2_rts.append(rt)
                    ms2_ints.append(sum(ints))

            all_rts = ms1_rts + ms2_rts
            ax[0].plot(ms1_rts, ms1_ints, label = 'g%s, s%s' %(groupi, samplei), marker = 'o')
            ax[1].plot(ms2_rts, ms2_ints, label = 'g%s, s%s' %(groupi, samplei), marker = 'o')

    ax[1].set_xlabel('Retention Time (min)')
    ax[0].set_ylabel('Intensity')
    ax[1].set_ylabel('Intensity')
    ax[0].set_title('MS1 TIC')
    ax[1].set_title('MS2 TIC')

    ax[0].ticklabel_format(style='sci',scilimits=(0,0),axis='y')
    ax[1].ticklabel_format(style='sci',scilimits=(0,0),axis='y')

    fig.savefig(os.path.join(options.out_dir, '%s_tic.jpg' %options.output_label), dpi = 300, bbox_inches = 'tight')
    return

def plot_acquisition_schema(options, run_template):

    fig = Figure()
    ax = fig.subplots(1,1)

    plot_time = 0
    for scani, scan in enumerate(run_template):
        if scan['order'] == 1:
            ax.add_patch(
                Rectangle(
                    (options.ms1_min_mz, plot_time),
                    options.ms1_max_mz - options.ms1_min_mz,
                    scan['length'], alpha = 0.5, facecolor = 'blue'
                )
            )

        if scan['order'] == 2:
            ax.add_patch(
                Rectangle(
                    (scan['isolation_range'][0], plot_time),
                    scan['isolation_range'][1] - scan['isolation_range'][0],
                    scan['length'], alpha = 0.5, facecolor = 'red'
                )
            )

        plot_time += scan['length']

    ax.legend(handles=[
        Patch(color = 'blue', label = 'MS1', alpha = 0.5),
        Patch(color = 'red', label = 'MS2', alpha = 0.5)
    ])

    ax.set_xlabel('m/z')
    ax.set_ylabel('Time (s)')
    ax.set_title('Acquisition Schema')
    ax.set_xlim((options.ms1_min_mz,options.ms1_max_mz))
    ax.set_ylim((0,plot_time))
    fig.savefig(os.path.join(options.out_dir, '%s_acquisition_schema.jpg' %options.output_label), dpi = 300, bbox_inches = 'tight')

    return

def plot_preview_graphics(options):

    fig = Figure()
    ax = fig.subplots(1,1)

    data = []
    for groupi in range(options.n_groups):
        for samplei in range(options.samples_per_group):

            mzml_file = os.path.join(
                options.out_dir, '%s_group_%s_sample_%s.mzML' %(
                    options.output_label, groupi, samplei
                )
            )

            for (rt, lvl, mzs, ints) in MZMLReader(mzml_file):
                if lvl != 1: continue
                data.extend([
                    {'rt': rt, 'mz': mzs[i], 'intensity': ints[i]} for i in range(len(ints)) if ints[i] > 0
                ])

    df = pd.DataFrame(data)

    # heatmap plot
    ax.scatter(x = df['rt'], y = df['mz'], c = df['intensity'])
    ax.set_xlabel('Retention Time (min)')
    ax.set_ylabel('m/z')
    ax.ticklabel_format(useOffset=False)
    ax.set_title('MS1 m/z vs retention time heatmap')
    fig.savefig(os.path.join(options.out_dir, 'preview_heatmap.jpg'), dpi = 300, bbox_inches = 'tight')

    # find scan with highest max intensity
    max_int_rt = df[df['intensity'] == df['intensity'].max()].iloc[0].rt
    scan_df = df[df['rt'] == max_int_rt]

    # plot max intensity spectrum
    fig = Figure()
    ax = fig.subplots(1,1)
    ax.plot(scan_df['mz'], scan_df['intensity'], marker = 'o')
    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')
    ax.set_title('MS1 mass spectrum')
    fig.savefig(os.path.join(options.out_dir, 'preview_ms1_scan.jpg'), dpi = 300, bbox_inches = 'tight')

    return
