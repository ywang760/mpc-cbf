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


