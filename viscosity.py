import argparse
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import curve_fit
from datetime import datetime


def fit_viscosity(acf, temperature, volume, fraction):
    """Fit viscosity according to alghoritm in [Zhang, 2015] paper.

    Parameters:
        acf : 3d ndarray. Autocorrelation functions values [P] 
              with corresponding times [ps] for every simulation.
        temperature : 1d ndarray. Temperature [K] values for every
                      simulation.
        volume : 1d ndarrray. Volume [nm^3] values for every 
                 simulation.
        fraction : threshold value for STD to average ratio.
    Returns:
        viscosity : Fitted value of viscosity [P].
        stats : Some statistics about fitting.
    """
    viscosities, time = get_viscosity_samples(acf, temperature, volume)
    viscosity_ave = np.average(viscosities, axis=0)
    viscosity_std = np.std(viscosities, axis=0)
    threshold_index = get_t_cut_index(viscosity_ave, viscosity_std, fraction)
    [A_std, b], std_cov = curve_fit(power_function, time, viscosity_std)
    time_trunc = time[0:threshold_index]
    viscosity_ave_trunc = viscosity_ave[0:threshold_index]
    weight = 1 / (time_trunc**b)
    initial_guess = [0, 0.5, 1, 1]
    bounds = ([0, 0, 0, 0], [np.inf, 1, np.inf, np.inf])
    [A, alpha, tau_1, tau_2], ave_cov = curve_fit(double_exponential_function, 
                                                  time_trunc, 
                                                  viscosity_ave_trunc,
                                                  sigma=weight,
                                                  absolute_sigma=True,
                                                  p0=initial_guess,
                                                  bounds=bounds,
                                                  max_nfev=1000
                                                 )
    viscosity = A * alpha * tau_1 + A * (1 - alpha) * tau_2
    stats = {"samples": viscosities.shape[0],
             "time": time,
             "average": viscosity_ave,
             "std": viscosity_std,
             "threshold_index": threshold_index,
             "t_cut": time[threshold_index],
             "A_std": A_std,
             "b": b,
             "std_cov": std_cov,
             "A": A,
             "alpha": alpha,
             "tau_1": tau_1,
             "tau_2": tau_2,
             "ave_cov": ave_cov,
             "fit": double_exponential_function(time, A, alpha, tau_1, tau_2)
            }
    return viscosity, stats
    

def power_function(t, A, b):
    "Return A*t^b."
    y = A * t**b
    return y


def double_exponential_function(t, A, alpha, tau_1, tau_2):
    y = A * (alpha * tau_1 * (1 - np.exp(-t / tau_1)) 
             + (1 - alpha) * tau_2 * (1 - np.exp(-t / tau_2)))
    return y


def get_t_cut_index(average, std, fraction):
    """Calculate ``t_cut`` index. 

    Parameters:
        average : 1d ndarray. Average values.
        std : 1d ndarray. STD values.
        fraction : threshold value for STD to average ratio.
    Returns:
        Index of first STD to average ratio value greater than 
        fraction.
    """
    std_vs_ave = std / average
    lower_than_fraction = np.argwhere(std_vs_ave >= fraction)
    if len(lower_than_fraction) > 0:
        index = lower_than_fraction[0, 0]
    else:
        index = std_vs_ave.shape[0] - 1
    return index


def load_data(project_directory, stem, filename, amount=0):
    """Load data from files in specified directories.
    Looks for files with ``filename`` in directories with common name
    part ``stem`` inside the ``project_directory`` and loads its data.

    Parameters:
        project_directory : Name of the directory with the 
                            simulations directories.
        stem : Common part of the simulations directories names.
        filename : Name of file with data.
        amount : Maximal number of the simulations directories for 
                 looking. If less or equal zero or greater than number
                 of all available directories look and load all 
                 possible data. Default is 0.
    Returns:
        ndarray with data. 
    """
    if not Path(project_directory).is_dir():
        raise OSError(project_directory + " doesn't exist.")
    data = []
    directories = list(Path(project_directory).glob(stem + "*"))
    if not directories:
        raise OSError("There is no " + stem 
                      + " in the " + project_directory + ".")
    if (amount > 0) and (amount < len(directories)):
        directories = directories[0:amount]
    for d in directories:
        path_to_file = d / Path(filename)
        data.append(np.loadtxt(path_to_file, comments=("#", "@", "&")))
    data = np.array(data)
    return data 


def get_viscosity_samples(acf, temperature, volume):
    """Obtain viscosity samples from autocorrelation data.

    Parameters:
        acf : 3d ndarray. Autocorrelation functions values [P] 
              with corresponding times [ps] for every simulation.
        temperature : 1d ndarray. Temperature [K] values for every
                      simulation.
        volume : 1d ndarrray. Volume [nm^3] values for every 
                 simulation.
    Returns:
        viscosities : 2d ndarray. Viscosity time dependencies values 
                      for every simulation [P].
        time : 1d ndarray. Time values for calculated viscosities [ps].
    """
    viscosities = []
    for a, t, v in zip(acf, temperature, volume): 
        viscosity = use_green_kubo(a, t[0], v) 
        viscosities.append(viscosity[:, 1])
    time = viscosity[:, 0]
    viscosities = np.array(viscosities)
    return viscosities, time


def use_green_kubo(acf, temperature, volume):
    """Calculate time dependence of viscosity from 
    Green-Kubo relation.

    Parameters:
        acf : 2xN ndarray. Time values [ps] and corresponding 
              autocorrelation function data [bar^2].
        temperature : Float. Temperature of the system [K]. 
        volume : Float. Volume of the simulation box [nm^3].
    Returns:
        2x(N-1) ndarray. Time values [p] and corresponding viscosity 
        time dependence values [P].
    """
    k_Boltzmann = 1.380649e-23
    unit_conversion_coefficient = 1e-28
    running_integral = cumulative_trapezoid(acf[:, 1], acf[:, 0])
    viscosity = (unit_conversion_coefficient * running_integral 
                 * volume / (k_Boltzmann * temperature))
    out = np.column_stack((acf[1:, 0], viscosity))
    return out


if __name__ == "__main__":
    start_time = datetime.now()
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", 
                        help=("The name of project directory. Default is "
                              "current directory"), 
                        nargs="?",
                        default=".")
    parser.add_argument("stem", 
                        help=("The common name part for directories "
                              "with the results of calculations. " 
                              "Default is traj_"),
                        nargs="?",
                        default="traj_"
                        )
    parser.add_argument("--fraction", 
                        help=("Threshold value for STD to average ratio. "
                              "Used for calculation of t_cut."
                              "Default is 0.4"),
                        default=0.4,
                        type=float
                        )
    parser.add_argument("--amount", 
                        help=("Number of directories used for analysis. "
                              "If lesser than 0 or greater than number "
                              "of available directories, all directories "
                              "will be used. Default is 0"),
                        default=0,
                        type=int
                        )
    parser.add_argument("--oa", 
                        help=("File name for saving the averaged viscosity "
                              "time dependence [ps, P]. If not set this data "
                              "will not be saved."), 
                        metavar="NAME",
                        )
    parser.add_argument("--os", 
                        help=("File name for saving the viscosity STD time "
                              "dependence [ps, P]. If not set this data "
                              "will not be saved."),
                        metavar="NAME",
                        )
    parser.add_argument("--of", 
                        help=("File name for saving the fitted viscosity "
                              "time dependence [ps, P]. If not set this "
                              "data will not be saved."),
                        metavar="NAME",
                        )
    args = parser.parse_args()
    acf = load_data(args.directory, args.stem, "acf.xvg", args.amount)
    temperature = load_data(args.directory, args.stem, "temperature.dat", 
                            args.amount)
    volume = load_data(args.directory, args.stem, "volume.dat", args.amount)
    viscosity_fit, stats = fit_viscosity(acf, 
                                         temperature, 
                                         volume, 
                                         args.fraction)
    if args.amount == 1:
        sample_addition_effect = "---"
    else:
        viscosity_fit_prev, stats_prev = fit_viscosity(acf[:-1], 
                                                       temperature[:-1], 
                                                       volume[:-1], 
                                                       args.fraction)
        sample_addition_effect = np.abs(viscosity_fit_prev - viscosity_fit)
    if args.oa:
        np.savetxt(args.oa, np.column_stack((stats["time"], stats["average"])))
    if args.os:
        np.savetxt(args.os, np.column_stack((stats["time"], stats["std"])))
    if args.of:
        np.savetxt(args.of, np.column_stack((stats["time"], stats["fit"])))
    run_time = datetime.now() - start_time
    format_string = "{0:35s} : {1}"
    print(format_string.format("Run time", run_time))
    print(format_string.format("Number of samples", stats["samples"]))
    print(format_string.format("t_cut, [ps]", stats["t_cut"]))
    print(format_string.format("t_cut index", stats["threshold_index"]))
    print(format_string.format("STD A", stats["A_std"]))
    print(format_string.format("STD b", stats["b"]))
    print(format_string.format("STD fitting condition number", 
                                 np.linalg.cond(stats["std_cov"])))
    print(format_string.format("STD fitting covariances", 
                                 np.diag(stats["std_cov"])))
    print(format_string.format("A", stats["A"]))
    print(format_string.format("alpha", stats["alpha"]))
    print(format_string.format("tau_1", stats["tau_1"]))
    print(format_string.format("tau_2", stats["tau_2"]))
    print(format_string.format("Average fitting condition number", 
                                 np.linalg.cond(stats["ave_cov"])))
    print(format_string.format("Average fitting covariances", 
                                 np.diag(stats["ave_cov"])))
    print(format_string.format("Sample addition effect", 
                                 sample_addition_effect))
    print(format_string.format("Viscosity, [P]", viscosity_fit))
