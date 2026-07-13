import numpy as np
import matplotlib.pyplot as plt
import corner
import matplotlib.lines as mlines

def overlapping_corner_three(array1, array2, array3, param_labels=None,sample_labels=None, save_plot=False, filename="corner_plot.pdf"):
    """
    Creates a corner plot using three 2D arrays, overlapping them with different colors,
    and displays histograms with consistent binning for comparison. Histograms are unfilled,
    the y-axis is kept on the first histogram, the x-axis is kept on the last histogram, and
    the last histogram has a label beneath the x-axis showing the parameter it represents.

    Parameters
    ----------
    array1, array2, array3 : 2D arrays
        The input datasets. Each should be of shape (n_samples, n_params), where `n_samples`
        is the number of samples, and `n_params` is the number of parameters. Arrays can
        have different lengths (n_samples).
    param_labels : list of str, optional
        The labels for each parameter to use in the corner plot. Defaults to None, in which
        case generic labels will be used.
    save_plot : bool, optional
        Whether to save the plot to a file. Defaults to False.
    filename : str, optional
        The filename to save the plot if `save_plot=True`. Defaults to "corner_plot.pdf".
    
    Returns
    -------
    None
        Displays the corner plot.
    """

    # Set up colors for the three datasets
    colors = ['red', 'blue', 'green']

    # Compute consistent bin edges based on combined data
    def get_consistent_bins(arrays, n_params, n_bins=50):
        bin_edges = []
        for i in range(n_params):
            combined_data = np.concatenate([array[:, i] for array in arrays])
            bins = np.histogram_bin_edges(combined_data, bins=n_bins)
            bin_edges.append(bins)
        return bin_edges

    # Compute bin edges for all parameters
    arrays = [array1, array2, array3]
    n_params = array1.shape[1]
    bin_edges = get_consistent_bins(arrays, n_params)

    # Create the corner plot for the first dataset (without 1D histograms)
    figure = corner.corner(array1, smooth=0.9,
                           label_kwargs=dict(fontsize=17),
                           title_kwargs=dict(fontsize=16),
                           levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
                           labels=param_labels, color=colors[0], show_titles=False, 
                           plot_datapoints=False, plot_density=True, fill_contours=True)
    
    # Overlap the second and third datasets onto the existing figure
    corner.corner(array2, smooth=0.9,
                  label_kwargs=dict(fontsize=17),
                  title_kwargs=dict(fontsize=16),
                  levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
                  labels=param_labels, color=colors[1], fig=figure, show_titles=False, 
                  plot_datapoints=False, plot_density=True, fill_contours=True)
    corner.corner(array3, smooth=0.9,
                  label_kwargs=dict(fontsize=17),
                  title_kwargs=dict(fontsize=16),
                  levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
                  labels=param_labels, color=colors[2], fig=figure, show_titles=False, 
                  plot_datapoints=False, plot_density=True, fill_contours=True)

    # Add histograms with consistent binning to the 1D marginal plots (diagonal panels)
    axes = np.array(figure.axes).reshape((n_params, n_params))
    for i in range(n_params):
        ax = axes[i, i]
        
        # Clear the existing diagonal histogram (which was automatically created by `corner`)
        ax.clear()
        param_name = param_labels[i] if param_labels else f"Param {i+1}"
        ax.set_title(param_name, fontsize=20, pad=75)  # Increase padding to make space for median/std
        
        # Plot unfilled histograms with consistent bins for each dataset on the diagonal
        for array, label, color, bins in zip(arrays, sample_labels, colors, [bin_edges[i]]*3):
            ax.hist(array[:, i], bins=bins, histtype='step', color=color, label=label, density=True)

        # Compute and add median and 68% quantile range in LaTeX format
        for idx, (array, color) in enumerate(zip(arrays, colors), start=1):
            median = np.percentile(array[:, i], 50)
            lower_quantile = np.percentile(array[:, i], 16)
            upper_quantile = np.percentile(array[:, i], 84)
            lower_diff = median - lower_quantile
            upper_diff = upper_quantile - median
            # Use LaTeX-like format for displaying the uncertainties
            ax.text(0.5, 0.9 + 0.17 * idx, 
                    f"${median:.3f}_{{-{lower_diff:.3f}}}^{{+{upper_diff:.3f}}}$", 
                    color=color, transform=ax.transAxes, fontsize=16, ha='center')

        # Keep y-axis visible for the first histogram, x-axis visible for the last histogram
        
        ax.set_yticks([])  # Hide y-axis for all other diagonal histograms

        if i == n_params - 1:
            ax.xaxis.set_visible(True)  # Keep x-axis visible for the last histogram
            ax.set_xlabel(param_labels[i]) # Add label beneath the x-axis for the last histogram
            plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
        else:
            ax.set_xticks([])  # Hide x-axis for all other diagonal histograms
    # Adjust layout: Reduce the spacing between the panels
    figure.subplots_adjust(top=0.9, wspace=0.05, hspace=0.05)

    # Add the legend
    handles = [
        mlines.Line2D([], [], color=colors[i], label=sample_labels[i]) for i in range(len(sample_labels))
    ]
    figure.legend(handles=handles, loc='upper right', bbox_to_anchor=(0.85, 0.85), fontsize=20, frameon=False)
    
    # Set tick label font size for non-diagonal plots
    label_font_size = 13
    for ax in figure.get_axes():
        ax.tick_params(axis='both', labelsize=label_font_size)
    ax = axes[-1, -1]
    ax.annotate(param_labels[-1], xy=(0.5, -0.38), xycoords='axes fraction', fontsize=20, ha='center', va='center')
    # Show or save the plot
    if save_plot:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Corner plot saved as {filename}.")
    else:
        plt.show()
        


def overlapping_corner_two(array1, array2, param_labels=None, sample_labels=None, save_plot=False, filename="corner_plot.pdf"):
    """
    Creates a corner plot using three 2D arrays, overlapping them with different colors,
    and displays histograms with consistent binning for comparison. Histograms are unfilled,
    the y-axis is kept on the first histogram, the x-axis is kept on the last histogram, and
    the last histogram has a label beneath the x-axis showing the parameter it represents.

    Parameters
    ----------
    array1, array2 : 2D arrays
        The input datasets. Each should be of shape (n_samples, n_params), where `n_samples`
        is the number of samples, and `n_params` is the number of parameters. Arrays can
        have different lengths (n_samples).
    param_labels : list of str, optional
        The labels for each parameter to use in the corner plot. Defaults to None, in which
        case generic labels will be used.
    save_plot : bool, optional
        Whether to save the plot to a file. Defaults to False.
    filename : str, optional
        The filename to save the plot if `save_plot=True`. Defaults to "corner_plot.pdf".
    
    Returns
    -------
    None
        Displays the corner plot.
    """

    # Set up colors for the three datasets
    colors = ['red', 'blue']

    # Compute consistent bin edges based on combined data
    def get_consistent_bins(arrays, n_params, n_bins=50):
        bin_edges = []
        for i in range(n_params):
            combined_data = np.concatenate([array[:, i] for array in arrays])
            bins = np.histogram_bin_edges(combined_data, bins=n_bins)
            bin_edges.append(bins)
        return bin_edges

    # Compute bin edges for all parameters
    arrays = [array1, array2,]
    n_params = array1.shape[1]
    bin_edges = get_consistent_bins(arrays, n_params)

    # Create the corner plot for the first dataset (without 1D histograms)
    figure = corner.corner(array1, smooth=0.9,
                           label_kwargs=dict(fontsize=17),
                           title_kwargs=dict(fontsize=16),
                           levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
                           labels=param_labels, color=colors[0], show_titles=False, 
                           plot_datapoints=False, plot_density=True, fill_contours=True)
    
    # Overlap the second and third datasets onto the existing figure
    corner.corner(array2, smooth=0.9,
                  label_kwargs=dict(fontsize=17),
                  title_kwargs=dict(fontsize=16),
                  levels=(1 - np.exp(-0.5), 1 - np.exp(-2), 1 - np.exp(-9 / 2.)),
                  labels=param_labels, color=colors[1], fig=figure, show_titles=False, 
                  plot_datapoints=False, plot_density=True, fill_contours=True)


    # Add histograms with consistent binning to the 1D marginal plots (diagonal panels)
    axes = np.array(figure.axes).reshape((n_params, n_params))
    for i in range(n_params):
        ax = axes[i, i]
        
        # Clear the existing diagonal histogram (which was automatically created by `corner`)
        ax.clear()
        param_name = param_labels[i] if param_labels else f"Param {i+1}"
        ax.set_title(param_name, fontsize=20, pad=55)  # Increase padding to make space for median/std
        
        # Plot unfilled histograms with consistent bins for each dataset on the diagonal
        for array, label, color, bins in zip(arrays, sample_labels, colors, [bin_edges[i]]*3):
            ax.hist(array[:, i], bins=bins, histtype='step', color=color, label=label, density=True)

        # Compute and add median and 68% quantile range in LaTeX format
        for idx, (array, color) in enumerate(zip(arrays, colors), start=1):
            median = np.percentile(array[:, i], 50)
            lower_quantile = np.percentile(array[:, i], 16)
            upper_quantile = np.percentile(array[:, i], 84)
            lower_diff = median - lower_quantile
            upper_diff = upper_quantile - median
            # Use LaTeX-like format for displaying the uncertainties
            ax.text(0.5, 0.9 + 0.17 * idx, 
                    f"${median:.3f}_{{-{lower_diff:.3f}}}^{{+{upper_diff:.3f}}}$", 
                    color=color, transform=ax.transAxes, fontsize=16, ha='center')

        # Keep y-axis visible for the first histogram, x-axis visible for the last histogram
        
        ax.set_yticks([])  # Hide y-axis for all other diagonal histograms

        if i == n_params - 1:
            ax.xaxis.set_visible(True)  # Keep x-axis visible for the last histogram
            ax.set_xlabel(param_labels[i]) # Add label beneath the x-axis for the last histogram
            plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
        else:
            ax.set_xticks([])  # Hide x-axis for all other diagonal histograms
        
    # Adjust layout: Reduce the spacing between the panels
    figure.subplots_adjust(top=0.9, wspace=0.05, hspace=0.05)

    # Add the legend
    handles = [
        mlines.Line2D([], [], color=colors[i], label=sample_labels[i]) for i in range(len(sample_labels))
    ]
    figure.legend(handles=handles, loc='upper right', bbox_to_anchor=(0.95, 0.95), fontsize=20, frameon=False)
    
    # Set tick label font size for non-diagonal plots
    label_font_size = 13
    for ax in figure.get_axes():
        ax.tick_params(axis='both', labelsize=label_font_size)
    ax = axes[-1, -1]
    ax.annotate(param_labels[-1], xy=(0.5, -0.38), xycoords='axes fraction', fontsize=20, ha='center', va='center')
    # Show or save the plot
    if save_plot:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Corner plot saved as {filename}.")
    else:
        plt.show()
