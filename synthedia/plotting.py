import os

from matplotlib.figure import Figure
from matplotlib.patches import Rectangle, Patch

from .mzml import MZMLReader
import numpy as np
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
            ax[0].plot(ms1_rts, ms1_ints, label = 'g%s, s%s' %(groupi, samplei))
            ax[1].plot(ms2_rts, ms2_ints, label = 'g%s, s%s' %(groupi, samplei))

    ax[1].set_xlabel('Retention Time (min)')
    ax[0].set_ylabel('Intensity')
    ax[1].set_ylabel('Intensity')
    ax[0].set_title('MS1')
    ax[1].set_title('MS2')
    fig.savefig(os.path.join(options.out_dir, '%s_tic.jpg' %options.output_label), dpi = 300, bbox_inches = 'tight')
    return

def plot_acquisition_schema(options, run_template):

    fig = Figure()
    ax = fig.subplots(1,1)

    plot_time = 0
    for i in range(3):
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


