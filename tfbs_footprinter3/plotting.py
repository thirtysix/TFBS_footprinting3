"""Promoter-level matplotlib plotting for tfbs_footprinter3.

Renders the final promoter figure (predicted TFBSs, conservation, CpG,
eQTLs, metaclusters, ATAC-Seq, CAGE, and TF-expression correlation tracks)
to SVG via matplotlib.

Matplotlib is a large (~1.5s) import; we defer it to the first call to
plot_promoter/plot_promoter_all via `_ensure_matplotlib_loaded()`. This
keeps `tfbs_footprinter3 --help`, `-no`, and HPC runs (which pass `-no`)
from paying the import cost.
"""
from __future__ import annotations

import math
import os
from operator import itemgetter

# Module-level placeholders populated by _ensure_matplotlib_loaded().
# Type as Any so callers aren't bothered by the None-init pattern.
plt = None  # matplotlib.pyplot
mpl = None  # pylab.mpl
numpyrandom = None  # numpy.random


def _ensure_matplotlib_loaded():
    """Lazily import matplotlib + pylab + numpy.random on first figure call.

    Matplotlib's Agg backend is selected BEFORE the pyplot import, so any
    prior import of pyplot in the process wins over this. In practice the
    tool never imports pyplot elsewhere.
    """
    global plt, mpl, numpyrandom
    if plt is not None:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as _plt
    from numpy import random as _numpyrandom
    from pylab import mpl as _mpl
    plt = _plt
    mpl = _mpl
    numpyrandom = _numpyrandom


def plot_promoter(target_species, transcript_id, species_group, alignment, alignment_len, promoter_before_tss, promoter_after_tss, transcript_name, top_x_greatest_hits_dict, target_dir, converted_reg_dict, converted_gerps_in_promoter, cpg_list, converted_cages, converted_metaclusters_in_promoter, converted_atac_seqs_in_promoter, converted_eqtls):
    """
    Plot the predicted TFBSs, onto a 5000 nt promoter graph, which possess support above the current strand threshold.
    ['binding_prot', 'species', 'motif', 'strand', 'start', 'end', 'TSS-relative start', 'TSS-relative end', 'PWM score', 'p-value', 'pos in align.', 'combined affinity score', 'support']
    """
    _ensure_matplotlib_loaded()

    # set axes for human
    if target_species == "homo_sapiens":
        fig = plt.figure(figsize=(10, 6))
        ax1 = plt.subplot2grid((20,1),(0,0), rowspan = 6, colspan = 11)
        ax8 = plt.subplot2grid((20,1),(6,0), rowspan = 2, colspan = 11)
        ax2 = plt.subplot2grid((20,1),(8,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax3 = plt.subplot2grid((20,1),(10,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax4 = plt.subplot2grid((20,1),(12,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax5 = plt.subplot2grid((20,1),(14,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax6 = plt.subplot2grid((20,1),(16,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax7 = plt.subplot2grid((20,1),(18,0), sharex=ax1, rowspan = 2, colspan = 11)

        # Set format of the plot(s)
        # Hide x-ticks for all plots except the lowest
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax3.get_xticklabels(), visible=False)
        plt.setp(ax4.get_xticklabels(), visible=False)
        plt.setp(ax5.get_xticklabels(), visible=False)
        plt.setp(ax6.get_xticklabels(), visible=False)
        plt.setp(ax8.get_xticklabels(), visible=False)

        # plt + ax labels
        ax1.text(1.02,.5,'Predicted TFBSs', verticalalignment='center', transform=ax1.transAxes, rotation='vertical', fontsize=8)
        ax1.set_ylabel("Combined Affinity Score", fontsize = 8, labelpad = 0)
        ax1.text(1.005,0.99,'+ strand', verticalalignment='top', transform=ax1.transAxes, rotation='vertical', fontsize=6)
        ax1.text(1.005,.01,'- strand', verticalalignment='bottom', transform=ax1.transAxes, rotation='vertical', fontsize=6)
        ax2.text(1.01,.5,'GERP\nConserv.\n'+species_group, verticalalignment='center', transform=ax2.transAxes, rotation='vertical', fontsize=5)
        ax3.text(1.01,.5,'CpG\nObs/Exp', verticalalignment='center', transform=ax3.transAxes, rotation='vertical', fontsize=6)
        ax4.text(1.01,.5,'eQTLs', verticalalignment='center', transform=ax4.transAxes, rotation='vertical', fontsize=6)
        ax5.text(1.01,.5,'TFBS\nMeta\nClusters', verticalalignment='center', transform=ax5.transAxes, rotation='vertical', fontsize=6)
        ax6.text(1.01,.5,'ATAC-Seq', verticalalignment='center', transform=ax6.transAxes, rotation='vertical', fontsize=6)
        ax7.text(1.01,.5,'CAGE\nPeaks\n(TSSs)', verticalalignment='center', transform=ax7.transAxes, rotation='vertical', fontsize=6)
        ax8.text(1.01,.5,'TF\nExpress.\nCorr.', verticalalignment='center', transform=ax8.transAxes, rotation='vertical', fontsize=6)

    ### as of now the data for non-human species is limited to predicted TFBSs, conservation, and CpG
    else:
        fig = plt.figure(figsize=(10, 6))
        ax1 = plt.subplot2grid((10,1),(0,0), rowspan = 6, colspan = 11)
        ax2 = plt.subplot2grid((10,1),(6,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax3 = plt.subplot2grid((10,1),(8,0), sharex=ax1, rowspan = 2, colspan = 11)

        # Set format of the plot(s)
        # Hide x-ticks for all plots except the lowest
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)

        # plt + ax labels
        ax1.text(1.02,.5,'Predicted TFBSs', verticalalignment='center', transform=ax1.transAxes, rotation='vertical', fontsize=8)
        ax1.set_ylabel("Combined Affinity Score", fontsize = 8, labelpad = 0)
        ax1.text(1.005,0.99,'+ strand', verticalalignment='top', transform=ax1.transAxes, rotation='vertical', fontsize=6)
        ax1.text(1.005,.01,'- strand', verticalalignment='bottom', transform=ax1.transAxes, rotation='vertical', fontsize=6)
        ax2.text(1.01,.5,'GERP\nConserv.\n'+species_group, verticalalignment='center', transform=ax2.transAxes, rotation='vertical', fontsize=5)
        ax3.text(1.01,.5,'CpG\nObs/Exp', verticalalignment='center', transform=ax3.transAxes, rotation='vertical', fontsize=6)


    # plot title
    title_str = target_species+"\n"+" ".join([transcript_name, transcript_id])
##    fig.text(0.065, 0.5, title_str, horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, rotation='vertical', fontsize=14)
    fig.text(0.065, 0.5, title_str, horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14)
    mpl.rcParams['axes.linewidth'] = 1.1
    plt.xlabel("Nucleotide position relative to TSS", labelpad=5)

    # Generate data for each of the greatest_hits and plot corresponding bar
    color_series=['#FFB300','#803E75','#FF6800','#A6BDD7','#C10020','#CEA262','#817066','#007D34','#F6768E','#00538A','#FF7A5C','#53377A','#FF8E00','#B32851','#F4C800','#7F180D','#93AA00','#593315','#F13A13','#232C16']
    color_dict = {'CTCF':'#FF0000', 'TBP':'#FF00FF'}
    y_range = []
    labels_used = []

    # Create a list sorted descending by combined affinity score, so that lower support hits that overlap can be seen.
    sorted_by_ca_list = []
    for TF, great_hits in top_x_greatest_hits_dict.items():
        for great_hit in great_hits:
            sorted_by_ca_list.append(great_hit)

    # ref-point
    sorted_by_ca_list = sorted(sorted_by_ca_list, key=itemgetter(10)) # sort by combined affinity p-val

    ### AX1: Predicted TFBSs
    for sorted_great_hit in sorted_by_ca_list:

        tf_name = sorted_great_hit[0]

        # choose a unique color for each tf_name
        if tf_name not in color_dict:
            pick = numpyrandom.randint(0, len(color_series) - 1)
            picked_color = color_series[pick]
            color_series.remove(picked_color)
            color_dict[tf_name] = picked_color
        else:
            picked_color = color_dict[tf_name]

        # if the label has been used, set label to "", otherwise labels will repeat in legend
        if tf_name in labels_used:
            lab = ""
        else:
            lab = tf_name
            labels_used.append(tf_name)

        edge_color = picked_color

        x_series = []
        y_series = []

        # ref-point
        #'binding prot.', 'motif', 'strand', 'start', 'end', 'TSS-relative start', 'TSS-relative end', 'frame score', 'p-value', 'combined\naffinity\nscore', 'species\nweights\nsum', 'cage\nweights\nsum', 'eqtls\nweights\nsum', 'atac\nweights\nsum', 'metacluster\nweights\nsum', 'cpg\nweight', 'corr.\nweight\nsum']
        binding_site_start = sorted_great_hit[5]
        binding_site_end = sorted_great_hit[6]
        combined_affinity = sorted_great_hit[9]
        binding_strand = int(sorted_great_hit[2])

        TF_center_point = float(binding_site_start + binding_site_end)/2
        TF_width = abs(binding_site_start - binding_site_end)
        x_series.append(TF_center_point)
        y_series.append(combined_affinity * binding_strand)
        y_range.append(combined_affinity)
        ax1.bar(x_series, y_series, facecolor = picked_color, edgecolor = edge_color, linewidth=1, alpha=0.9, align = 'center', width = TF_width, label = lab)

    # Set y-axis height based on number of entries in alignment
    y_range.sort()
    tens_y = int(int(y_range[-1])/10) + 1

    # Ensembl regulatory information
    # All will be horizontally plotted in some shade of red
    if len(converted_reg_dict) > 0:
        alpha_gradient = 1.0
        alpha_gradient_dict = {1:0}
        for i in range(2,100):
            alpha_gradient_dict[i] = 1./i
        reg_height = 1
        reg_height = (tens_y)/4
        for reg_id, data in converted_reg_dict.items():
            converted_start = int(data['converted_start'])
            converted_end = int(data['converted_end'])
            # limit length to first two words so legend isn't overrun
            description = data['description']
            reg_x_series = []
            reg_y_series = []
            center_point = float(converted_start + converted_end)/2
            reg_x_series.append(center_point)
            reg_y_series.append(reg_height)
            reg_x_series.append(center_point)
            reg_y_series.append(reg_height * -1)
            reg_width = abs(converted_start - converted_end)
            ax1.bar(reg_x_series, reg_y_series, facecolor='red', edgecolor='red', alpha=alpha_gradient, align = 'center', width=reg_width, label=description)
            alpha_gradient -= alpha_gradient_dict[len(converted_reg_dict)]
            reg_height += 0.5

    ax1.axhline(0, color = 'black', linewidth=0.5)

    ### AX2: Add GERP conservation bars
    for converted_gerp_in_promoter in converted_gerps_in_promoter:
        converted_gerp_start = converted_gerp_in_promoter[0]
        converted_gerp_end = converted_gerp_in_promoter[1]
        alpha_gradient = 1
        gerp_height = 1

        gerp_x_series = []
        gerp_y_series = []
        gerp_midpoint = float(converted_gerp_start + converted_gerp_end)/2
        gerp_x_series.append(gerp_midpoint)
        gerp_y_series.append(gerp_height)

        gerp_width = abs(converted_gerp_start - converted_gerp_end)
        ax2.bar(gerp_x_series, gerp_y_series, facecolor='black', edgecolor='black', alpha=alpha_gradient, align = 'center', width=gerp_width)

    ax2.set_yticks([0, 1])
    plt.setp(ax2.get_yticklabels(), fontsize=6)
    ax2.set_ylim(0, 1)

    ### AX3: CpG plot
    # [1 C, 1 if G, 1 if CPG, CorG, num_cpg, obs2exp]
    obs2exp = [x[5] for x in cpg_list]
    ax3.plot(range(-1 * alignment_len + promoter_after_tss, promoter_after_tss), obs2exp, color = 'red')
    gpc = []
    top_obs2exp = ax3.get_ylim()[-1]

    for x in cpg_list:
        if x[2] == 0:
            gpc.append(x[2])
        else:
            if top_obs2exp <= 1:
                gpc.append(1)
            else:
                gpc.append(top_obs2exp)

    ax3.bar(range(-1 * alignment_len + promoter_after_tss, promoter_after_tss), gpc, color = 'black')
    if top_obs2exp <1:
        top_obs2exp = 1
    ax3.set_ylim(0, top_obs2exp)
    ax3.set_yticks([0, 0.6, 1])
    ax3.set_yticklabels([0, 0.6, 1], va='center')
    plt.setp(ax3.get_yticklabels(), fontsize=6)
    ax3.axhline(0.6, color = 'black', alpha = 0.4)

    ### human-specific experimental data
    if target_species == "homo_sapiens":

        ### AX7: CAGE plot
        cage_height = 1
        cage_labels = []
        for converted_cage in converted_cages:
            converted_cage_start = converted_cage[0]
            converted_cage_end = converted_cage[1]
            description = converted_cage[2]

            cage_x_series = []
            cage_y_series = []
            cage_center_point = float(converted_cage_start + converted_cage_end)/2
            cage_x_series.append(cage_center_point)
            cage_y_series.append(cage_height)

            cage_width = abs(converted_cage_start - converted_cage_end)
            ax7.bar(cage_x_series, cage_y_series, facecolor='black', edgecolor='black', align = 'center', width=cage_width, label=description)

            # add label for the CAGE peak
            if -1 * promoter_before_tss <= converted_cage_start <= promoter_after_tss + 1 or -1 * promoter_before_tss <= converted_cage_end <= promoter_after_tss + 1:
                plt.text(cage_center_point, cage_height, description, color="red", rotation = 270, fontsize=5, horizontalalignment='center', verticalalignment='top')

        ax7.axes.get_yaxis().set_visible(False)

        ### AX5: GTRD plot
        gtrd_height = 1
        for converted_metacluster_in_promoter in converted_metaclusters_in_promoter:
            converted_metacluster_start = converted_metacluster_in_promoter[0]
            converted_metacluster_end = converted_metacluster_in_promoter[1]
            metacluster_peak_count = converted_metacluster_in_promoter[2]
            alpha_gradient = 0.5 + (metacluster_peak_count/1220.0)/2

            gtrd_x_series = []
            gtrd_y_series = []
            gtrd_center_point = float(converted_metacluster_start + converted_metacluster_end)/2
            gtrd_x_series.append(gtrd_center_point)
            gtrd_y_series.append(gtrd_height)

            gtrd_width = abs(converted_metacluster_start - converted_metacluster_end)
            ax5.bar(gtrd_x_series, gtrd_y_series, facecolor='black', edgecolor='black', alpha=alpha_gradient, align = 'center', width=gtrd_width)
        ax5.axes.get_yaxis().set_visible(False)

        ### AX6: ATAC-Seq plot
        for converted_atac_seq_in_promoter in converted_atac_seqs_in_promoter:
            converted_atac_seq_start = converted_atac_seq_in_promoter[0]
            converted_atac_seq_end = converted_atac_seq_in_promoter[1]
            atac_seq_peak_score = converted_atac_seq_in_promoter[2]
            alpha_gradient = 0.5 + atac_seq_peak_score/93.234864

            gtrd_x_series = []
            gtrd_y_series = []
            gtrd_midpoint = float(converted_atac_seq_start + converted_atac_seq_end)/2
            gtrd_x_series.append(gtrd_midpoint)
            gtrd_y_series.append(gtrd_height)

            gtrd_width = abs(converted_atac_seq_start - converted_atac_seq_end)
            ax6.bar(gtrd_x_series, gtrd_y_series, facecolor='black', edgecolor='black', alpha=alpha_gradient, align = 'center', width=gtrd_width)
        ax6.axes.get_yaxis().set_visible(False)

        ### AX4: eQTLs plot
        colors = ["green", "red"]
        magnitudes = []
        for converted_eqtl in converted_eqtls:
            converted_eqtl_start, converted_eqtl_end, converted_eqtl_mag = converted_eqtl
            if -1 * promoter_before_tss <= converted_eqtl_start <= promoter_after_tss + 1 or -1 * promoter_before_tss <= converted_eqtl_end <= promoter_after_tss + 1:
                eqtl_midpoint = float(converted_eqtl_start + converted_eqtl_end)/2
                eqtl_width = abs(converted_eqtl_start - converted_eqtl_end)
                eqtl_x_series = []
                eqtl_y_series = []
                eqtl_x_series.append(eqtl_midpoint)
                eqtl_y_series.append(converted_eqtl_mag)
                magnitudes.append(converted_eqtl_mag)
                if converted_eqtl_mag > 0:
                    c = colors[0]
                else:
                    c = colors[1]
                ax4.bar(eqtl_x_series, eqtl_y_series, facecolor=c, edgecolor=c, align = 'center', width=eqtl_width)
    ##            # arrow does not format properly, perhaps due to size.  y value starts not at 0, and arrow wraps over itself.
    ##            ax4.arrow(eqtl_midpoint, 0, 0, converted_eqtl_mag, color=c, length_includes_head = True, lw=10, width=0.01)

        ax4_yticks = [-1,0,1]
        if len(magnitudes) > 0:
            magnitudes.sort()
            ax4_yticks = [math.floor(magnitudes[0]), 0, math.ceil(magnitudes[-1])]
        ax4.set_yticks(ax4_yticks)
        ax4.set_yticklabels(ax4_yticks, va='center')
        plt.setp(ax4.get_yticklabels(), fontsize=6)
        ax4.axhline(0.0, color = 'black', alpha = 0.4)

        ### AX8: cage_correlations
        # rebuild dict to have just the top correlation
        # ref-point
        plot_tfs_corrs_colors = [(tf_name, hits_list[0][16],  color_dict[tf_name]) for tf_name, hits_list in top_x_greatest_hits_dict.items()]
        plot_tfs_corrs_colors_sorted = sorted(plot_tfs_corrs_colors, key=itemgetter(1), reverse=True)
        ax8.bar(range(0, len(plot_tfs_corrs_colors_sorted)), [x[1] for x in plot_tfs_corrs_colors_sorted], color=[x[2] for x in plot_tfs_corrs_colors_sorted], edgecolor = "none")
        ax8.set_ylim(0, plot_tfs_corrs_colors_sorted[0][1]+1)
        ax8.set_xlim(-1, len(top_x_greatest_hits_dict))
        ax8.set_yticks([0, math.ceil(plot_tfs_corrs_colors_sorted[0][1])+1])
        plt.setp(ax8.get_yticklabels(), fontsize=6)


    ## set ticks
    # based on 100's
    if y_range[-1] <= 100:
        for falling_y_thresh in range(100, -1, -10):
            if y_range[-1] < falling_y_thresh:
                y_thresh = falling_y_thresh
        ax1.set_yticks(range(-1* y_thresh, y_thresh+1, 10))
    else:
        ax1.set_yticks(range(-1 * (((tens_y*10)/100)+1)*100, (((tens_y*10)/100)+2)*100, 100))

    # format y-labels
    ylabs=ax1.get_yticks().tolist()
    ylabs=[abs(x) for x in ylabs]
    ax1.set_yticklabels(ylabs)
    plt.setp(ax1.get_yticklabels(), fontsize=8)

    # Misc
    plt.xlim([-1 * promoter_before_tss, promoter_after_tss + 1])

    # legend
    num_cols = 6
    ax1.legend(bbox_to_anchor=[0., 1.1, 1.0, .102], loc='center', ncol=num_cols, prop={'size':8}, mode="expand", borderaxespad=0.)

    # produce .svg figure
    plt.subplots_adjust(hspace=0.40)
    fig.savefig(os.path.join(target_dir, os.path.basename(target_dir) + '.Promoterhisto'  + '.svg'), facecolor='white', bbox_inches='tight')
    plt.clf()
    plt.close()


def plot_promoter_all(target_species, transcript_id, species_group, alignment, alignment_len, promoter_before_tss, promoter_after_tss, transcript_name, top_x_greatest_hits_dict, target_dir, converted_reg_dict, converted_gerps_in_promoter, cpg_list, converted_cages, converted_metaclusters_in_promoter, converted_atac_seqs_in_promoter, converted_eqtls):
    """
    Plot the predicted TFBSs, onto a 5000 nt promoter graph, which possess support above the current strand threshold.
    ['binding_prot', 'species', 'motif', 'strand', 'start', 'end', 'TSS-relative start', 'TSS-relative end', 'PWM score', 'p-value', 'pos in align.', 'combined affinity score', 'support']
    """
    _ensure_matplotlib_loaded()

    if target_species == "homo_sapiens":
        fig = plt.figure(figsize=(10, 6))
        ax1 = plt.subplot2grid((20,1),(0,0), rowspan = 6, colspan = 11)
        ax8 = plt.subplot2grid((20,1),(6,0), rowspan = 2, colspan = 11)
        ax2 = plt.subplot2grid((20,1),(8,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax3 = plt.subplot2grid((20,1),(10,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax4 = plt.subplot2grid((20,1),(12,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax5 = plt.subplot2grid((20,1),(14,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax6 = plt.subplot2grid((20,1),(16,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax7 = plt.subplot2grid((20,1),(18,0), sharex=ax1, rowspan = 2, colspan = 11)

        # Set format of the plot(s)
        # Hide x-ticks for all plots except the lowest
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax3.get_xticklabels(), visible=False)
        plt.setp(ax4.get_xticklabels(), visible=False)
        plt.setp(ax5.get_xticklabels(), visible=False)
        plt.setp(ax6.get_xticklabels(), visible=False)
        plt.setp(ax8.get_xticklabels(), visible=False)

        # plt + ax labels
        ax1.text(1.02,.5,'Predicted TFBSs', verticalalignment='center', transform=ax1.transAxes, rotation='vertical', fontsize=8)
        ax1.set_ylabel("Combined Affinity Score", fontsize = 8, labelpad = 0)
        ax1.text(1.005,0.99,'+ strand', verticalalignment='top', transform=ax1.transAxes, rotation='vertical', fontsize=6)
        ax1.text(1.005,.01,'- strand', verticalalignment='bottom', transform=ax1.transAxes, rotation='vertical', fontsize=6)
        ax2.text(1.01,.5,'GERP\nConserv.\n'+species_group, verticalalignment='center', transform=ax2.transAxes, rotation='vertical', fontsize=5)
        ax3.text(1.01,.5,'CpG\nObs/Exp', verticalalignment='center', transform=ax3.transAxes, rotation='vertical', fontsize=6)
        ax4.text(1.01,.5,'eQTLs', verticalalignment='center', transform=ax4.transAxes, rotation='vertical', fontsize=6)
        ax5.text(1.01,.5,'TFBS\nMeta\nClusters', verticalalignment='center', transform=ax5.transAxes, rotation='vertical', fontsize=6)
        ax6.text(1.01,.5,'ATAC-Seq', verticalalignment='center', transform=ax6.transAxes, rotation='vertical', fontsize=6)
        ax7.text(1.01,.5,'CAGE\nPeaks\n(TSSs)', verticalalignment='center', transform=ax7.transAxes, rotation='vertical', fontsize=6)
        ax8.text(1.01,.5,'TF\nExpress.\nCorr.', verticalalignment='center', transform=ax8.transAxes, rotation='vertical', fontsize=6)

    ### as of now the data for non-human species is limited to predicted TFBSs, conservation, and CpG
    else:
        fig = plt.figure(figsize=(10, 6))
        ax1 = plt.subplot2grid((10,1),(0,0), rowspan = 6, colspan = 11)
        ax2 = plt.subplot2grid((10,1),(6,0), sharex=ax1, rowspan = 2, colspan = 11)
        ax3 = plt.subplot2grid((10,1),(8,0), sharex=ax1, rowspan = 2, colspan = 11)

        # Set format of the plot(s)
        # Hide x-ticks for all plots except the lowest
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)

        # plt + ax labels
        ax1.text(1.02,.5,'Predicted TFBSs', verticalalignment='center', transform=ax1.transAxes, rotation='vertical', fontsize=8)
        ax1.set_ylabel("Combined Affinity Score", fontsize = 8, labelpad = 0)
        ax1.text(1.005,0.99,'+ strand', verticalalignment='top', transform=ax1.transAxes, rotation='vertical', fontsize=6)
        ax1.text(1.005,.01,'- strand', verticalalignment='bottom', transform=ax1.transAxes, rotation='vertical', fontsize=6)
        ax2.text(1.01,.5,'GERP\nConserv.\n'+species_group, verticalalignment='center', transform=ax2.transAxes, rotation='vertical', fontsize=5)
        ax3.text(1.01,.5,'CpG\nObs/Exp', verticalalignment='center', transform=ax3.transAxes, rotation='vertical', fontsize=6)


    # plot title
    title_str = target_species+"\n"+" ".join([transcript_name, transcript_id])
    fig.text(0.065, 0.5, title_str, horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, rotation='vertical', fontsize=14)
    mpl.rcParams['axes.linewidth'] = 1.1
    plt.xlabel("Nucleotide position relative to TSS", labelpad=5)

    # Generate data for each of the greatest_hits and plot corresponding bar
    color_series=['#FFB300','#803E75','#FF6800','#A6BDD7','#C10020','#CEA262','#817066','#007D34','#F6768E','#00538A','#FF7A5C','#53377A','#FF8E00','#B32851','#F4C800','#7F180D','#93AA00','#593315','#F13A13','#232C16']
    color_dict = {'CTCF':'#FF0000', 'TBP':'#FF00FF'}
    y_range = []
    labels_used = []

    # Create a list sorted descending by combined affinity score, so that lower support hits that overlap can be seen.
    sorted_by_ca_list = []
    for TF, great_hits in top_x_greatest_hits_dict.items():
        for great_hit in great_hits:
            sorted_by_ca_list.append(great_hit)
    # ref-point
    sorted_by_ca_list = sorted(sorted_by_ca_list, key=itemgetter(9), reverse=True)

    ### AX1: Predicted TFBSs
    for sorted_great_hit in sorted_by_ca_list:
        tf_name = sorted_great_hit[0]

##        # choose a unique color for each tf_name
##        if tf_name not in color_dict:
##            pick = numpyrandom.randint(0, len(color_series) - 1)
##            picked_color = color_series[pick]
##            color_series.remove(picked_color)
##            color_dict[tf_name] = picked_color
##        else:
##            picked_color = color_dict[tf_name]

        picked_color = '#FFB300'

        # if the label has been used, set label to "", otherwise labels will repeat in legend
        if tf_name in labels_used:
            lab = ""
        else:
            lab = tf_name
            labels_used.append(tf_name)

        edge_color = picked_color

        x_series = []
        y_series = []

        # ref-point
        #'binding prot.', 'motif', 'strand', 'start', 'end', 'TSS-relative start', 'TSS-relative end', 'frame score', 'p-value', 'combined\naffinity\nscore', 'species\nweights\nsum', 'cage\nweights\nsum', 'eqtls\nweights\nsum', 'atac\nweights\nsum', 'metacluster\nweights\nsum', 'cpg\nweight', 'corr.\nweight\nsum']
        binding_site_start = sorted_great_hit[5]
        binding_site_end = sorted_great_hit[6]
        combined_affinity = sorted_great_hit[9]
        binding_strand = int(sorted_great_hit[2])

        TF_center_point = float(binding_site_start + binding_site_end)/2
        TF_width = abs(binding_site_start - binding_site_end)
        x_series.append(TF_center_point)
        y_series.append(combined_affinity * binding_strand)
        y_range.append(combined_affinity)
        ax1.bar(x_series, y_series, facecolor = picked_color, edgecolor = edge_color, linewidth=1, alpha=0.9, align = 'center', width = TF_width, label = lab)

    # Set y-axis height based on number of entries in alignment
    y_range.sort()
    tens_y = int(y_range[-1])/10 + 1

    # Ensembl regulatory information
    # All will be horizontally plotted in some shade of red
    if len(converted_reg_dict) > 0:
        alpha_gradient = 1.0
        alpha_gradient_dict = {1:0}
        for i in range(2,100):
            alpha_gradient_dict[i] = 1./i
        reg_height = 1
        reg_height = (tens_y)/4
        for reg_id, data in converted_reg_dict.items():
            converted_start = int(data['converted_start'])
            converted_end = int(data['converted_end'])
            # limit length to first two words so legend isn't overrun
            description = data['description']
            reg_x_series = []
            reg_y_series = []
            center_point = float(converted_start + converted_end)/2
            reg_x_series.append(center_point)
            reg_y_series.append(reg_height)
            reg_x_series.append(center_point)
            reg_y_series.append(reg_height * -1)
            reg_width = abs(converted_start - converted_end)
            ax1.bar(reg_x_series, reg_y_series, facecolor='red', edgecolor='red', alpha=alpha_gradient, align = 'center', width=reg_width, label=description)
            alpha_gradient -= alpha_gradient_dict[len(converted_reg_dict)]
            reg_height += 0.5


    ax1.axhline(0, color = 'black', linewidth=0.5)

    ### AX2: Add GERP conservation bars
    for converted_gerp_in_promoter in converted_gerps_in_promoter:
        converted_gerp_start = converted_gerp_in_promoter[0]
        converted_gerp_end = converted_gerp_in_promoter[1]
        alpha_gradient = 1
        gerp_height = 1

        gerp_x_series = []
        gerp_y_series = []
        gerp_midpoint = float(converted_gerp_start + converted_gerp_end)/2
        gerp_x_series.append(gerp_midpoint)
        gerp_y_series.append(gerp_height)

        gerp_width = abs(converted_gerp_start - converted_gerp_end)
        ax2.bar(gerp_x_series, gerp_y_series, facecolor='black', edgecolor='black', alpha=alpha_gradient, align = 'center', width=gerp_width)

    ax2.set_yticks([0, 1])
    plt.setp(ax2.get_yticklabels(), fontsize=6)
    ax2.set_ylim(0, 1)


    ### AX3: CpG plot
    # [1 C, 1 if G, 1 if CPG, CorG, num_cpg, obs2exp]
    obs2exp = [x[5] for x in cpg_list]
    ax3.plot(range(-1 * alignment_len + promoter_after_tss, promoter_after_tss), obs2exp, color = 'red')
    gpc = []
    top_obs2exp = ax3.get_ylim()[-1]

    for x in cpg_list:
        if x[2] == 0:
            gpc.append(x[2])
        else:
            if top_obs2exp <= 1:
                gpc.append(1)
            else:
                gpc.append(top_obs2exp)

    ax3.bar(range(-1 * alignment_len + promoter_after_tss, promoter_after_tss), gpc, color = 'black')
    if top_obs2exp <1:
        top_obs2exp = 1
    ax3.set_ylim(0, top_obs2exp)
    ax3.set_yticks([0, 0.6, 1])
    ax3.set_yticklabels([0, 0.6, 1], va='center')
    plt.setp(ax3.get_yticklabels(), fontsize=6)
    ax3.axhline(0.6, color = 'black', alpha = 0.4)

    ### human-specific experimental data
    if target_species == "homo_sapiens":

        ### AX7: CAGE plot
        cage_height = 1
        cage_labels = []
        for converted_cage in converted_cages:
            converted_cage_start = converted_cage[0]
            converted_cage_end = converted_cage[1]
            description = converted_cage[2]
##            if ".." in description:
##                description = ""
            cage_x_series = []
            cage_y_series = []
            cage_center_point = float(converted_cage_start + converted_cage_end)/2
            cage_x_series.append(cage_center_point)
            cage_y_series.append(cage_height)

            cage_width = abs(converted_cage_start - converted_cage_end)
            ax7.bar(cage_x_series, cage_y_series, facecolor='black', edgecolor='black', align = 'center', width=cage_width, label=description)

            # add label for the CAGE peak
            if -1 * promoter_before_tss <= converted_cage_start <= promoter_after_tss + 1 or -1 * promoter_before_tss <= converted_cage_end <= promoter_after_tss + 1:
                plt.text(cage_center_point, cage_height, description, color="red", rotation = 270, fontsize=5, horizontalalignment='center', verticalalignment='top')

        ax7.axes.get_yaxis().set_visible(False)

        ### AX5: GTRD plot
        gtrd_height = 1
        for converted_metacluster_in_promoter in converted_metaclusters_in_promoter:
            converted_metacluster_start = converted_metacluster_in_promoter[0]
            converted_metacluster_end = converted_metacluster_in_promoter[1]
            metacluster_peak_count = converted_metacluster_in_promoter[2]
            alpha_gradient = 0.5 + (metacluster_peak_count/1220.0)/2

            gtrd_x_series = []
            gtrd_y_series = []
            gtrd_center_point = float(converted_metacluster_start + converted_metacluster_end)/2
            gtrd_x_series.append(gtrd_center_point)
            gtrd_y_series.append(gtrd_height)

            gtrd_width = abs(converted_metacluster_start - converted_metacluster_end)
            ax5.bar(gtrd_x_series, gtrd_y_series, facecolor='black', edgecolor='black', alpha=alpha_gradient, align = 'center', width=gtrd_width)
        ax5.axes.get_yaxis().set_visible(False)

        ### AX6: ATAC-Seq plot
        for converted_atac_seq_in_promoter in converted_atac_seqs_in_promoter:
            converted_atac_seq_start = converted_atac_seq_in_promoter[0]
            converted_atac_seq_end = converted_atac_seq_in_promoter[1]
            atac_seq_peak_score = converted_atac_seq_in_promoter[2]
            alpha_gradient = 0.5 + atac_seq_peak_score/93.234864

            gtrd_x_series = []
            gtrd_y_series = []
            gtrd_midpoint = float(converted_atac_seq_start + converted_atac_seq_end)/2
            gtrd_x_series.append(gtrd_midpoint)
            gtrd_y_series.append(gtrd_height)

            gtrd_width = abs(converted_atac_seq_start - converted_atac_seq_end)
            ax6.bar(gtrd_x_series, gtrd_y_series, facecolor='black', edgecolor='black', alpha=alpha_gradient, align = 'center', width=gtrd_width)
        ax6.axes.get_yaxis().set_visible(False)

        ### AX4: eQTLs plot
        colors = ["green", "red"]
        magnitudes = []
        for converted_eqtl in converted_eqtls:
            converted_eqtl_start, converted_eqtl_end, converted_eqtl_mag = converted_eqtl
            if -1 * promoter_before_tss <= converted_eqtl_start <= promoter_after_tss + 1 or -1 * promoter_before_tss <= converted_eqtl_end <= promoter_after_tss + 1:
                eqtl_midpoint = float(converted_eqtl_start + converted_eqtl_end)/2
                eqtl_width = abs(converted_eqtl_start - converted_eqtl_end)
                eqtl_x_series = []
                eqtl_y_series = []
                eqtl_x_series.append(eqtl_midpoint)
                eqtl_y_series.append(converted_eqtl_mag)
                magnitudes.append(converted_eqtl_mag)
                if converted_eqtl_mag > 0:
                    c = colors[0]
                else:
                    c = colors[1]
                ax4.bar(eqtl_x_series, eqtl_y_series, facecolor=c, edgecolor=c, align = 'center', width=eqtl_width)
    ##            # arrow does not format properly, perhaps due to size.  y value starts not at 0, and arrow wraps over itself.
    ##            ax4.arrow(eqtl_midpoint, 0, 0, converted_eqtl_mag, color=c, length_includes_head = True, lw=10, width=0.01)

        ax4_yticks = [-1,0,1]
        if len(magnitudes) > 0:
            magnitudes.sort()
            ax4_yticks = [math.floor(magnitudes[0]), 0, math.ceil(magnitudes[-1])]
        ax4.set_yticks(ax4_yticks)
        ax4.set_yticklabels(ax4_yticks, va='center')
        plt.setp(ax4.get_yticklabels(), fontsize=6)
        ax4.axhline(0.0, color = 'black', alpha = 0.4)

    ## set ticks
    # based on 100's
    if y_range[-1] <= 100:
        for falling_y_thresh in range(100, -1, -10):
            if y_range[-1] < falling_y_thresh:
                y_thresh = falling_y_thresh
        ax1.set_yticks(range(-1* y_thresh, y_thresh+1, 10))
    else:
        ax1.set_yticks(range(-1 * (((tens_y*10)/100)+1)*100, (((tens_y*10)/100)+2)*100, 100))

    ylabs=ax1.get_yticks().tolist()
    ylabs=[abs(x) for x in ylabs]
    ax1.set_yticklabels(ylabs)
    plt.setp(ax1.get_yticklabels(), fontsize=8)

    # Misc
    plt.xlim([-1 * promoter_before_tss, promoter_after_tss + 1])

    # legend
    num_cols = 6

    # produce .svg figure
    plt.subplots_adjust(hspace=0.40)
    fig.savefig(os.path.join(target_dir, os.path.basename(target_dir) + '.Promoterhisto'  + '.all.svg'), facecolor='white', bbox_inches='tight')
    plt.clf()
    plt.close()

