import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator

def vertical_result_comparison(
    x_min, x_max, observable, results, theory_bands, output_name,
    group_separators=None, layout_mode='right', second_theory=False,
        output_extension='pdf',width=10,height=None,
):
    """
    Draw vertical comparison of results and save to file.

    Args:
        x_min (float): Minimum x-axis value.
        x_max (float): Maximum x-axis value.
        observable (str): Observable name / x-axis label.
        results (list): List of lists with format ["Name", value, uncertainty, color, marker_style, True].
        theory_bands (list): List of [label, value, uncertainty, color] for theory bands.
        output_name (str): Output filename (without extension).
        output_extension (str): File extension (default: pdf).
        group_separators (list or None): Indices after which to draw horizontal dashed lines.
        second_theory (bool): If True, adds label to the second listed theory band.
        layout_mode (str): One of 'right', 'below', 'outside', 'outside_below'.
    """
    outside_offset = (x_max - x_min) * 0.02
    label_offset = (x_max - x_min) * 0.01
    value_offset = (x_max - x_min) * 0.005
    n = len(results)
    if height is None:
        height=0.3*n+1
    fig, ax = plt.subplots(figsize=(width, height))

    ax.set_xlim(x_min, x_max)
    ax.set_xticks(np.linspace(x_min, x_max, 6)[1:-1])
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(axis='x', labelsize=14, direction="in")
    ax.tick_params(axis='x', which="minor", length=3, direction="in")
    #ax.xaxis.set_minor_locator(AutoMinorLocator())

    ax.set_ylim(-0.5, n - 0.5)
    ax.set_xlabel(observable, fontsize=14)

    ax.set_yticks([])
    ax.tick_params(axis='y', which='both', left=False)

    # Draw theory bands
    for i, (label, val, err, color) in enumerate(theory_bands):
        ax.axvspan(val - err, val + err, color=color, alpha=0.85, label=label if i == 0 else None)
        x_mid = (val - err + val + err) / 2
        if second_theory:
            if i == 1:
                ax.text(x_mid+.0001, n - 0.4, f"{val:.5f} ± {err:.5f} {label}", ha='right', va='bottom', fontsize=13, color=color)
            else:
                ax.text(x_mid-0.0001, n - 0.4, f"{label} {val:.5f} ± {err:.5f}", ha='left', va='bottom', fontsize=13, color=color)
        else:
            if i == 0:
                # Align label outside if requested
                if layout_mode in ['outside']:
                    ax.text(x_max + outside_offset, n - 0.4, f"{val:.5f} ± {err:.5f}",
                            ha='left', va='bottom', fontsize=13, color=color)
                else:
                    ax.text(x_max, n - 0.4, f"{val:.5f} ± {err:.5f}",
                            ha='right', va='bottom', fontsize=13, color=color)

    # Draw result points and texts
    for i, (name, value, uncertainty, color, marker_style, is_public) in enumerate(reversed(results)):
        y = i
        ax.errorbar(value, y, xerr=uncertainty, fmt=marker_style, markersize=7, color=color, linewidth=3)

        if layout_mode == 'right':
            ax.text(x_min + label_offset, y, name, ha='left', va='center', fontsize=13, fontweight='bold')
            ax.text(x_max - value_offset, y, f"{value:.5f} ± {uncertainty:.5f}", ha='right', va='center', fontsize=13)

        elif layout_mode == 'below':
            ax.text(x_min + label_offset, y + 0.1, name, ha='left', va='center', fontsize=13, fontweight='bold')
            ax.text(x_min + label_offset + (x_max - x_min) * 0.02, y - 0.15,
                    f"{value:.5f} ± {uncertainty:.5f}", ha='left', va='top', fontsize=13)

        elif layout_mode == 'outside':
            ax.text(x_min - (x_max - x_min) * 0.02, y, name,
                    ha='right', va='center', fontsize=14, fontweight='bold',color=color)
            ax.text(x_max + (x_max - x_min) * 0.02, y,
                    f"{value:.5f} ± {uncertainty:.5f}", ha='left', va='center', fontsize=14,color=color)

        elif layout_mode == 'outside_below':
            ax.text(x_min - (x_max - x_min) * 0.02, y + 0.1,
                    name, ha='right', va='center', fontsize=13, fontweight='bold')
            ax.text(x_min - (x_max - x_min) * 0.02, y - 0.2,
                    f"{value:.5f} ± {uncertainty:.5f}", ha='right', va='top', fontsize=13)

        else:
            raise ValueError(f"Unknown layout_mode: {layout_mode}. Use 'right', 'below', 'outside', or 'outside_below'.")

    # Draw horizontal dashed lines for group separators
    if group_separators:
        for sep_index in group_separators:
            y_sep = n - sep_index - 0.5
            ax.axhline(y=y_sep, linestyle='--', color='gray', linewidth=1, alpha=0.5)

    ax.set_title("")
    fig.tight_layout()
    plt.savefig(f"{output_name}.{output_extension}")
    plt.close()



if __name__ == "__main__":
    # Example usage based on results shown in SMP-22-010
    results_ee = [
        ["LEP $A_{FB}^{0,l}$",    0.23099, 0.00053, 'black',        'o', True],
        ["LEP $P_{\\tau}$",         0.23159, 0.00041, 'black',        'o', True],
        ["LEP $A_{FB}^{0,b}$", 0.23221, 0.00029, 'black',    'o', True],
        ["LEP $A_{FB}^{0,c}$", 0.23220, 0.00081, 'black',    'o', True],
        ["LEP $Q_{FB}^{had}$", 0.2324, 0.0012, 'black',    'o', True],
        ["SLD $A_{LR}^{0}$",          0.23098, 0.00026, 'black',        'o', True],
        ["LEP+SLD",          0.23153, 0.00016, 'black',        'o', True],
    ]

    results = [
        ["LEP $A_{FB}^{0,b}$", 0.23221, 0.00029, 'black',    'o', True],
        ["SLD $A_{LR}^{0}$",          0.23098, 0.00026, 'black',        'o', True],
        ["CDF 1.96 TeV",           0.23221, 0.00046, '#228B22',      'o', True],  # forest green
        ["D0 1.96 TeV",            0.23095, 0.00040, '#228B22',      'o', True],
        ["ATLAS 7 TeV",         0.23080, 0.00120, 'blue',         'o', True],
        ["LHCb 7+8 TeV",        0.23142, 0.00106, 'blue',         'o', True],
        ["CMS 8 TeV",           0.23101, 0.00053, 'blue',         'o', True],
        #["ATLAS 8 TeV",         0.23140, 0.00036, 'blue',         'o', False],  # isPreliminary
        ["LHCb 13 TeV",         0.23147, 0.00050, 'blue',         'o', True], 
        ["CMS 13 TeV",          0.23152, 0.00031, 'blue',          'o', True],
        #["This analysis 13 TeV",          0.23156, 0.00024, 'red',          'o', True],
        ["CMS 13 TeV + Aux.",          0.23156, 0.00024, 'red',          'o', True],
        ["CMS 13+13.6 TeV\n+ Aux. (expected)",          0.23161, 0.00020, 'red',          'o', True],
    ]

    theory_bands = [
        ["SM", 0.23161, 0.00004, '#FFA500'],  
    ]

    xmin=0.229;
    xmax=0.234;
    observable = r"$\sin^{2}\theta_{\mathrm{eff}}^{\ell}$"

    # vertical_result_comparison(xmin, xmax, observable, results, theory_bands, "my_test_v1", group_separators=[4, 6, 10], layout_mode='right')
    # vertical_result_comparison(0.228, 0.233, observable, results, theory_bands, "my_test_v2", group_separators=[4, 6, 10], layout_mode='below')
    # vertical_result_comparison(0.229, 0.233, observable, results, theory_bands, "my_test_v3", group_separators=[4, 6, 10], layout_mode='outside')
    # vertical_result_comparison(0.229, 0.233, observable, results, theory_bands, "my_test_v4", group_separators=[4, 6, 10], layout_mode='outside_below')

    # Example usage with second theory band
    theory_bands = [
        ["SM", 0.23161, 0.00004, '#FFA500'],
        ["2HDM", 0.23110, 0.00010, "#1EA5FF"  ]
    ]

    # vertical_result_comparison(xmin, xmax, observable, results, theory_bands, "my_test_v1", group_separators=[4, 6, 10], layout_mode='right', second_theory=True)
    # vertical_result_comparison(0.228, 0.233, observable, results, theory_bands, "my_test_v2", group_separators=[4, 6, 10], layout_mode='below', second_theory=True)
    # vertical_result_comparison(0.229, 0.233, observable, results, theory_bands, "my_test_v3", group_separators=[4, 6, 10], layout_mode='outside', second_theory=True)
    # vertical_result_comparison(0.229, 0.233, observable, results, theory_bands, "my_test_v4", group_separators=[4, 6, 10], layout_mode='outside_below', second_theory=True)

    ### s2eff paper
    #vertical_result_comparison(xmin, xmax, observable, results, theory_bands, "fig/Fig7_all_sw2", group_separators=[2, 4, 9], layout_mode='outside', second_theory=True,width=8,height=6)
    #vertical_result_comparison(xmin, xmax, observable, results_ee, theory_bands, "fig/Fig1_LEP_SLD", group_separators=[6,], layout_mode='outside', second_theory=True,width=8,height=4)

    #vertical_result_comparison(xmin, xmax, observable, results, theory_bands, "test", group_separators=[2, 4, 9], layout_mode='outside', second_theory=True,width=8,height=6)


    results_leptonphoton = [
        ["LEP $A_{FB}^{0,l}$",    0.23099, 0.00053, 'black',        'o', True],
        ["LEP $P_{\\tau}$",         0.23159, 0.00041, 'black',        'o', True],
        ["LEP $A_{FB}^{0,b}$", 0.23221, 0.00029, 'black',    'o', True],
        ["LEP $A_{FB}^{0,c}$", 0.23220, 0.00081, 'black',    'o', True],
        ["LEP $Q_{FB}^{had}$", 0.2324, 0.0012, 'black',    'o', True],
        ["SLD $A_{LR}^{0}$",          0.23098, 0.00026, 'black',        'o', True],
        ["CDF 1.96 TeV",           0.23221, 0.00046, '#228B22',      'o', True],  # forest green
        ["D0 1.96 TeV",            0.23095, 0.00040, '#228B22',      'o', True],
        ["ATLAS 7 TeV",         0.23080, 0.00120, 'blue',         'o', True],
        ["LHCb 7+8 TeV",        0.23142, 0.00106, 'blue',         'o', True],
        ["CMS 8 TeV",           0.23101, 0.00053, 'blue',         'o', True],
        ["LHCb 13 TeV",         0.23147, 0.00050, 'blue',         'o', True],
        ["CMS 13 TeV",          0.23152, 0.00031, 'blue',          'o', True],
        ["This analysis 13 TeV",          0.23156, 0.00024, 'red',          'o', True],
    ]
        
    results_ewprecision = [
        # ["LEP $A_{FB}^{0,l}$",    0.23099, 0.00053, 'black',        'o', True],
        # ["LEP $P_{\\tau}$",         0.23159, 0.00041, 'black',        'o', True],
        # ["LEP $A_{FB}^{0,b}$", 0.23221, 0.00029, 'black',    'o', True],
        # ["LEP $A_{FB}^{0,c}$", 0.23220, 0.00081, 'black',    'o', True],
        # ["LEP $Q_{FB}^{had}$", 0.2324, 0.0012, 'black',    'o', True],
        # ["SLD $A_{LR}^{0}$",          0.23098, 0.00026, 'black',        'o', True],
        # ["CDF 1.96 TeV",           0.23221, 0.00046, '#228B22',      'o', True],  # forest green
        # ["D0 1.96 TeV",            0.23095, 0.00040, '#228B22',      'o', True],
        ["ATLAS 7 TeV",         0.23080, 0.00120, 'blue',         'o', True],
        ["LHCb 7+8 TeV",        0.23142, 0.00106, 'blue',         'o', True],
        ["CMS 8 TeV",           0.23101, 0.00053, 'blue',         'o', True],
        ["ATLAS 8 TeV (preliminary)",         0.23140, 0.00036, 'blue',         'o', False],  # isPreliminary
        ["LHCb 13 TeV",         0.23147, 0.00050, 'blue',         'o', True], 
        ["CMS 13 TeV",          0.23152, 0.00031, 'blue',          'o', True],
        ["Improved CMS 13 TeV",          0.23156, 0.00024, '#228B22',          'o', True],
    ]

    theory_bands = [
        ["SM", 0.23161, 0.00004, '#FFA500'],  
    ]
    #vertical_result_comparison(xmin, xmax, observable, results_ewprecision, theory_bands, "sw2", group_separators=[6, 8, 15], layout_mode='outside', second_theory=True,width=10,height=6)
    #vertical_result_comparison(xmin, xmax, observable, results_ewprecision, theory_bands, "sw2", group_separators=[7], layout_mode='outside', second_theory=True,width=10,height=3)

    # For SMP-25-012 DY+b
    theory_bands = [
        ["SM", 0.23161, 0.00004, '#FFA500'],
    ]

    results_won = [
        ["LEP $A_{FB}^{0,l}$", 0.23099, 0.00053, 'black', 'o', True],
        ["LEP $P_{\\tau}$",    0.23159, 0.00041, 'black', 'o', True],
        ["LEP $A_{FB}^{0,b}$", 0.23221, 0.00029, 'black', 'o', True],
        ["LEP $A_{FB}^{0,c}$", 0.23220, 0.00081, 'black', 'o', True],
        ["LEP $Q_{FB}^{had}$", 0.2324,  0.0012,  'black', 'o', True],
        ["SLD $A_{LR}^{0}$",   0.23098, 0.00026, 'black', 'o', True],
        ["CDF 1.96 TeV",       0.23221, 0.00046, '#228B22', 'o', True],  # forest green
        ["D0 1.96 TeV",        0.23095, 0.00040, '#228B22', 'o', True],
        ["ATLAS 7 TeV",        0.2308,  0.0012,  'blue', 'o', True],
        ["LHCb 7+8 TeV",       0.23142, 0.00106, 'blue', 'o', True],
        ["CMS 8 TeV",          0.23101, 0.00053, 'blue', 'o', True],
        ["CMS 13 TeV",         0.23152, 0.00031, 'blue', 'o', True],
        ["This analysis 13 TeV", 0.23495, 0.00292, 'red', 'o', True],
    ]

    xmin = 0.229;
    xmax = 0.239;
    vertical_result_comparison(xmin, xmax, observable, results_won, theory_bands, "sin2w_table", group_separators=[6, 8, 11], layout_mode='outside', second_theory=True, width=8, height=6)

    results_won = [
        ["LEP $A_{FB}^{0,l}$", 0.23099, 0.00053, 'black', 'o', True],
        ["LEP $P_{\\tau}$",    0.23159, 0.00041, 'black', 'o', True],
        ["LEP $A_{FB}^{0,b}$", 0.23221, 0.00029, 'black', 'o', True],
        ["LEP $A_{FB}^{0,c}$", 0.23220, 0.00081, 'black', 'o', True],
        ["LEP $Q_{FB}^{had}$", 0.2324,  0.0012,  'black', 'o', True],
        ["SLD $A_{LR}^{0}$",   0.23098, 0.00026, 'black', 'o', True],
        ["CDF 1.96 TeV",       0.23221, 0.00046, '#228B22', 'o', True],  # forest green
        ["D0 1.96 TeV",        0.23095, 0.00040, '#228B22', 'o', True],
        ["ATLAS 7 TeV",        0.2308,  0.0012,  'blue', 'o', True],
        ["LHCb 7+8 TeV",       0.23142, 0.00106, 'blue', 'o', True],
        ["CMS 8 TeV",          0.23101, 0.00053, 'blue', 'o', True],
        ["CMS 13 TeV",         0.23152, 0.00031, 'blue', 'o', True],
        ["This analysis",      0.23211, 0.00263, 'red', 'o', True],
        ["ee",                 0.22930, 0.00433, 'red', 'o', True],
        ["$\mu\mu$",           0.23369, 0.00328, 'red', 'o', True],
    ]

    xmin = 0.224;
    xmax = 0.238;
    vertical_result_comparison(xmin, xmax, observable, results_won, theory_bands, "sin2w_table_backup_cos1p0", group_separators=[6, 8, 12], layout_mode='outside', second_theory=True, width=8, height=6)

    results_won = [
        ["LEP $A_{FB}^{0,l}$", 0.23099, 0.00053, 'black', 'o', True],
        ["LEP $P_{\\tau}$",    0.23159, 0.00041, 'black', 'o', True],
        ["LEP $A_{FB}^{0,b}$", 0.23221, 0.00029, 'black', 'o', True],
        ["LEP $A_{FB}^{0,c}$", 0.23220, 0.00081, 'black', 'o', True],
        ["LEP $Q_{FB}^{had}$", 0.2324,  0.0012,  'black', 'o', True],
        ["SLD $A_{LR}^{0}$",   0.23098, 0.00026, 'black', 'o', True],
        ["CDF 1.96 TeV",       0.23221, 0.00046, '#228B22', 'o', True],  # forest green
        ["D0 1.96 TeV",        0.23095, 0.00040, '#228B22', 'o', True],
        ["ATLAS 7 TeV",        0.2308,  0.0012,  'blue', 'o', True],
        ["LHCb 7+8 TeV",       0.23142, 0.00106, 'blue', 'o', True],
        ["CMS 8 TeV",          0.23101, 0.00053, 'blue', 'o', True],
        ["CMS 13 TeV",         0.23152, 0.00031, 'blue', 'o', True],
        ["This analysis",      0.23495, 0.00292, 'red', 'o', True],
        ["ee",                 0.23074, 0.00468, 'red', 'o', True],
        ["$\mu\mu$",           0.23754, 0.00364, 'red', 'o', True],
    ]

    xmin = 0.225;
    xmax = 0.242;
    vertical_result_comparison(xmin, xmax, observable, results_won, theory_bands, "sin2w_table_backup_cos0p9", group_separators=[6, 8, 12], layout_mode='outside', second_theory=True, width=8, height=6)
