import numpy as np
import matplotlib.pyplot as plt

def CI_compute(sample_arr):
    '''
    compute the 95% confidence interval for samples.
    :param sample_arr: [entries, M] steps is the number of statistics entry, M is the number of samples.
    :return:
    mean_arr [entries,]
    ci [entries,]
    '''
    # constants
    T, M = sample_arr.shape
    # compute the mean_arr
    mean_arr = np.mean(sample_arr, axis=1)  # [T,]
    # compute the std_arr
    std_arr = np.std(sample_arr, axis=1)  # [T,]
    # compute the 95% CI
    ci = 1.96 * (std_arr / np.sqrt(M))  # [T,]
    return mean_arr, ci

def percentile_compute(sample_arr):
    '''
    Compute the median, 25th percentile (Q1), and 75th percentile (Q3) for samples.

    :param sample_arr: [entries, M] where entries is the number of statistics entry,
                       M is the number of samples.
    :return:
    median_arr [entries,]
    q1_arr [entries,] (25th percentile)
    q3_arr [entries,] (75th percentile)
    '''
    # Compute the median
    median_arr = np.median(sample_arr, axis=1)  # [T,]

    # Compute the 25th percentile (Q1)
    q1_arr = np.percentile(sample_arr, 25, axis=1)  # [T,]
    q1_arr = median_arr - q1_arr

    # Compute the 75th percentile (Q3)
    q3_arr = np.percentile(sample_arr, 75, axis=1)  # [T,]
    q3_arr = q3_arr - median_arr

    return median_arr, q1_arr, q3_arr

def CI_compute_with_inf(sample_arr):
    # constants
    T, M = sample_arr.shape
    valid_sample_arr = [[] for i in range(T)]
    mean_arr = []
    ci = []
    for i in range(T):
        for j in range(M):
            if sample_arr[i, j] != float("inf"):
                valid_sample_arr[i].append(sample_arr[i, j])

    for i in range(T):
        valid_size = len(valid_sample_arr[i])
        if valid_size > 0:
            mean = np.mean(valid_sample_arr[i])
            mean_arr.append(mean)
            std = np.std(valid_sample_arr[i])  # [T,]
            c = 1.96 * (std / np.sqrt(valid_size))
            ci.append(c)
        else:
            mean_arr.append(np.nan)
            ci.append(np.nan)

    return np.array(mean_arr), np.array(ci)

def percentile_compute_with_inf(sample_arr):
    # constants
    T, M = sample_arr.shape
    valid_sample_arr = [[] for i in range(T)]
    median_arr = []
    q1_arr = []
    q3_arr = []
    for i in range(T):
        for j in range(M):
            if sample_arr[i, j] != float("inf"):
                valid_sample_arr[i].append(sample_arr[i, j])

    for i in range(T):
        valid_size = len(valid_sample_arr[i])
        if valid_size > 0:
            median = np.median(valid_sample_arr[i])
            median_arr.append(median)
            q1 = np.percentile(valid_sample_arr[i], 25)
            q1 = median - q1
            # Compute the 75th percentile (Q3)
            q3 = np.percentile(valid_sample_arr[i], 75)
            q3 = q3 - median
            q1_arr.append(q1)
            q3_arr.append(q3)
        else:
            median_arr.append(np.nan)
            q1_arr.append(np.nan)
            q3_arr.append(np.nan)

    return np.array(median_arr), np.array(q1_arr), np.array(q3_arr)

def CI_plot(x, mean_arr, ci, save_name="./ci_plot.png", xlabel="entry", ylabel="value", label = ""):
    '''
    plot the confident interval for the reward curve across the training epochs.
    :param
    x: [entries,],
    mean_arr: [entries,], entries of properties of curve.
    ci: [entries,]
    :return: confidence interval plot for value v.s. entries.
    '''
    assert x.shape == mean_arr.shape
    # constants
    T, = mean_arr.shape
    # create canvas
    fig, ax = plt.subplots()
    # plot
    ax.plot(x, mean_arr[:], label=label)
    ax.fill_between(x,
                    mean_arr[:] - ci[:],
                    mean_arr[:] + ci[:],
                    alpha=0.1)
    plt.xticks(x)

    ####################################################
    # # assistant plot
    # assistant_x = np.linspace(1, 10, 100)
    # assistant_y = 0.13/assistant_x
    # ax.plot(assistant_x, assistant_y, ls = "--")
    ####################################################


    # legend and label, save the figure
    plt.legend()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig(save_name, bbox_inches='tight')

def list_CI_plot(x, mean_arr_list, negative_ci_arr_list, positive_ci_arr_list, label_list, colors, xlabel="entry", ylabel="value", save_name="./ci_plot.png"):
    '''
    plot the confident interval for the reward curve across the training epochs.
    :param
    x: [entries,],
    mean_arr: [entries,], entries of properties of curve.
    ci: [entries,]
    :return: confidence interval plot for value v.s. entries.
    '''
    assert x.shape == mean_arr_list[0].shape
    # create canvas
    fig, ax = plt.subplots(figsize=(10, 3))

    linestyles = ["--", "--", "--", "-", "-", "-", "-.", "-."]
    for i in range(len(mean_arr_list)):
    # plot
        ax.plot(x, mean_arr_list[i][:], c=colors[i], label=f"${label_list[i]}$", ls=linestyles[i])
        ax.scatter(x, mean_arr_list[i][:], c=colors[i], s=4, lw=1)
        ax.fill_between(x,
                        mean_arr_list[i][:] - negative_ci_arr_list[i][:],
                        mean_arr_list[i][:] + positive_ci_arr_list[i][:],
                        color=colors[i] ,alpha=0.1)
    plt.xticks(x, fontsize=13)
    plt.yticks(fontsize=13)

    ####################################################
    # # assistant plot
    # assistant_x = np.linspace(1, 10, 100)
    # assistant_y = 0.13/assistant_x
    # ax.plot(assistant_x, assistant_y, ls = "--")
    ####################################################


    # legend and label, save the figure
    plt.legend(fontsize=9, ncol=3)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    plt.savefig(save_name, bbox_inches='tight')

def histogram_list_CI_plot(categories, samples_arr_list, mean_arr_list, negative_err_arr_list, positive_err_arr_list, label_list, colors, xlabel="entry", ylabel="value", legend=False, save_name="./hist_plot.png"):
    n_category = len(categories)
    n_subgroup = len(mean_arr_list)
    values = mean_arr_list
    negative_cis = negative_err_arr_list
    positive_cis = positive_err_arr_list
    samples = samples_arr_list

    # hatches = ["/","/","/", None, None, None, None, None]
    alphas = [0.3, 0.3, 0.3, 0.7, 0.7, 0.7, 0.9, 0.9]
    ls = ["--", "--", "--", "-", "-", "-", "-.", "-.",]

    # values = [[] for i in range(n_category)]
    # cis = [[] for i in range(n_category)]
    # samples = [[] for i in range(n_category)]
    # for i in range(n_category):
    #     for j in range(n_subgroup):
    #         values[i].append(mean_arr_list[j][i])
    #         cis[i].append(ci_arr_list[j][i])
    #         samples[i].append(samples_arr_list[j][i])

    subgroups = label_list

    # Parameters for the bar chart
    num_categories = len(categories)
    num_subgroups = len(subgroups)
    width = 0.09  # Width of each bar

    # create canvas
    fig, ax = plt.subplots(figsize=(10, 3))
    x = np.arange(len(categories)) + categories[0]
    for i in range(num_subgroups):
        # Bar plot
        bar_positions = x + i * width - width * (num_subgroups - 1) / 2
        ax.bar(
            bar_positions,
            values[i],
            width,
            yerr=np.array([negative_cis[i], positive_cis[i]]),
            capsize=1.5,
            label=f"${subgroups[i]}$",
            color=colors[i],
            edgecolor="k",
            alpha=alphas[i],
            # hatch=hatches[i],
            ls = ls[i],
            lw=0.8
        )
        # Add scatter points with hollow circles
        for j, pos in enumerate(bar_positions):
            ax.scatter(
                [pos] * len(samples[i][j]),
                samples[i][j],
                facecolors='none',  # Hollow center
                edgecolors='black',  # Black border
                alpha=0.6,
                s=4,  # Point size
                lw=0.5
            )

    # legend and label, save the figure
    # plt.legend()
    if legend:
        ax.legend(
            fontsize=12,
            bbox_to_anchor=(0.47, 1.15),  # Center it above the plot
            loc='lower center',          # Place it at the bottom of the legend box
            ncol=3,                      # Display the legend in a single row
            borderaxespad=-1              # Remove padding between the axes and the legend
        )

    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_xticks(categories)
    ax.set_xticklabels(categories, fontsize=13)
    plt.yticks(fontsize=13)
    ax.yaxis.grid(True, linestyle='--', alpha=0.6)
    plt.savefig(save_name, bbox_inches='tight')


def Multi_CI_plot(x_axis_list, mu_list, ci_list, metric_names, xlabel=None, save_name="./NumberScale/NumberScaleMultiPlot"):
    assert len(metric_names) == 2  # TODO change this number for different number of CI plots
    # plot the overlap for different exp
    fig, ax1 = plt.subplots()
    fig.set_size_inches(14, 10)

    ax1.set_xlabel(xlabel, fontsize=30)
    ax1.tick_params(axis='x', labelsize=28)
    ax1.set_ylabel(metric_names[0], color="blue", fontsize=30, labelpad=10)
    # ax1.set_ylim(40, 80)
    line1 = ax1.plot(x_axis_list, mu_list[0], label=metric_names[0], color="blue", linestyle='-', lw=5, marker="o", ms=15)
    ax1.fill_between(x_axis_list,
                     (mu_list[0] - ci_list[0]),
                     (mu_list[0] + ci_list[0]),
                     alpha=0.1, color="blue")
    ax1.tick_params(axis='y', labelcolor='blue', color="blue", labelsize=28)

    ax2 = ax1.twinx()
    ax2.set_ylabel(metric_names[1], color="green", fontsize=30, labelpad=10)
    # ax2.set_ylabel(metric_names[1], color="green", fontsize=30)
    line2 = ax2.plot(x_axis_list, mu_list[1], label=metric_names[1], color="green", linestyle='--', lw=5, marker="^", ms=15)
    ax2.fill_between(x_axis_list,
                     (mu_list[1] - ci_list[1]),
                     (mu_list[1] + ci_list[1]),
                     alpha=0.1, color="green")
    # ax2.set_ylim([0, 30])
    ax2.tick_params(axis='y', labelcolor='green', color="green", labelsize=28)

    # ##############################################################################
    # ax3 = ax1.twinx()
    # # Offset the right spine of twin2.  The ticks and label have already been
    # # placed on the right by twinx above.
    # ax3.spines.right.set_position(("axes", 1.1))
    # ax3.set_ylabel(metric_names[2], color="red", fontsize=30, labelpad=10)
    # # ax3.set_ylabel(metric_names[2], color="red", fontsize=30)
    # line3 = ax3.plot(x_axis_list, mu_list[2], label=metric_names[2], color="red", linestyle='-.', lw=5, marker="s", ms=15)
    # ax3.fill_between(x_axis_list,
    #                  (mu_list[2] - ci_list[2]),
    #                  (mu_list[2] + ci_list[2]),
    #                  alpha=0.1, color="red")
    # ax3.tick_params(axis='y', labelcolor='red', color="red", labelsize=28)
    # ##############################################################################

    # Solution for having two legends
    legend = line1 + line2
    labs = [l.get_label() for l in legend]
    # ax1.legend(legend, labs, loc="lower right", fontsize=25, ncol=3)
    ax1.legend(legend, labs, loc="upper right", fontsize=25)

    # fig.tight_layout()
    plt.xticks(x_axis_list)
    # plt.title("Metrics")
    plt.savefig(save_name, bbox_inches='tight')


