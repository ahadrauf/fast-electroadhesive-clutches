import numpy as np
import matplotlib.pyplot as plt
from import_strain_test_data import *
from scipy.signal import butter, sosfiltfilt
from scipy.stats import linregress


def calculate_Fpreload(file_name, w_pin, depth_pin, verbose=False):
    time_loadcell, load_cell, relay_input, relay_output, idx_stop_moving, idx_relay_input_rise, idx_relay_input_fall = import_data(file_name,
                                                                                                                                   zero_load_cell=False,
                                                                                                                                   convert_to_kPa=False,
                                                                                                                                   convert_time_to_seconds=False)
    mu_objet_s = 0.281  # 0.173  # 0.28
    mu_pvdf_s = 0.188  # 0.154  # 0.38
    mu_objet_k = 0.173
    mu_pvdf_k = 0.154
    mu_tot_k = mu_objet_k + mu_pvdf_k
    # mu_tot_s = mu_objet_s + mu_pvdf_s
    initial_friction = np.nanmean(load_cell[idx_relay_input_rise - int(12800*0.5):idx_relay_input_rise])
    m_pin = 8338.53879*w_pin*depth_pin*100e-3 + 0.725823081e-3  # from 20241106 - TMech - Mass of Pins.xlsx (g --> kg)
    m_pvdf = 1.78e3*24e-6*60e-3*w_pin  # https://piezopvdf.com/CTFE07-terpolymer-film/
    g = 9.81
    F_preload = (1/mu_tot_k)*(initial_friction - m_pin*g*mu_objet_k - m_pvdf*g*mu_tot_k)
    if verbose:
        print("Initial Friction:", initial_friction, m_pin*g*mu_objet_k, m_pvdf*g*mu_tot_k)
    return F_preload


def calculate_engagement_time(file_name, plot=False):  # this method works pretty well?
    time_loadcell, load_cell, relay_input, relay_output, idx_stop_moving, idx_relay_input_rise, idx_relay_input_fall = import_data(file_name,
                                                                                                                                   zero_load_cell=False,
                                                                                                                                   convert_to_kPa=False,
                                                                                                                                   convert_time_to_seconds=False)
    sos = butter(5, (48, 52), btype='bandstop', fs=12800, output='sos')
    load_cell = sosfiltfilt(sos, load_cell)
    sos = butter(3, 250, btype='lowpass', fs=12800, output='sos')
    load_cell = sosfiltfilt(sos, load_cell)

    engagement_times = []
    r2scores = []
    idx_relay_input_rise = max(np.where(relay_output > 50)[0][0], np.where(relay_input > 0.5)[0][0])

    # Check whether the load cell engaged weirdly early (this part has no effect on the rest of the function)
    offset_when_voltage_enabled = np.nanmean(load_cell[idx_relay_input_rise - int(12800*0.5):idx_relay_input_rise - int(12800*0.002)])
    std_when_voltage_enabled = np.nanstd(load_cell[idx_relay_input_rise - int(12800*0.5):idx_relay_input_rise - int(12800*0.002)])
    idx_lc_23 = idx_relay_input_rise*2//3
    try:
        if np.where(load_cell[idx_lc_23:] > offset_when_voltage_enabled + 10*std_when_voltage_enabled)[0][0] + idx_lc_23 < idx_relay_input_rise:
            idx_load_cell_rise = np.where(load_cell[idx_lc_23:] > offset_when_voltage_enabled + 3*std_when_voltage_enabled)[0][0] + idx_lc_23
            print("The load cell weirdly engaged at", time_loadcell[idx_load_cell_rise], "before voltage was enabled at",
                  time_loadcell[idx_relay_input_rise], "Mean:", offset_when_voltage_enabled, "Std:", std_when_voltage_enabled,
                  "Threshold:", offset_when_voltage_enabled + 10*std_when_voltage_enabled)
            return np.infty
    except IndexError as e:
        print("The load cell weirdly engaged, and threw an error in the error checking code, returning infinity")
        return np.infty

    offset_when_voltage_enabled = np.nanmean(load_cell[idx_relay_input_rise - int(12800*0.01):idx_relay_input_rise - int(12800*0.002)])
    didx_fit_start = int(12800*0.0005)
    didx_fits = np.arange(int(12800*0.002), int(12800*0.02))
    for didx_fit in didx_fits:
        time_fit = time_loadcell[idx_relay_input_rise + didx_fit_start:idx_relay_input_rise + didx_fit]
        slope, intercept, r, p, se = linregress(time_fit, load_cell[idx_relay_input_rise + didx_fit_start:idx_relay_input_rise + didx_fit])
        engagement_time = ((offset_when_voltage_enabled - intercept)/slope - time_loadcell[idx_relay_input_rise])
        engagement_times.append(engagement_time)
        r2scores.append(r*r)
    engagement_times = np.array(engagement_times)
    r2scores = np.array(r2scores)
    engagement_times = engagement_times[r2scores >= 0.8]
    # r2scores = r2scores[(r2scores >= 0.8)]
    engagement_times_filtered = engagement_times[engagement_times >= 0]
    print("% of Engagement Times after Filtering", len(engagement_times_filtered)/len(engagement_times))  # , "Cutoff", r2score_cutoff)
    if len(engagement_times_filtered) > 0:
        engagement_time = np.min(engagement_times_filtered)
    else:
        print("No engagement times > 0, Max Time:", np.max(engagement_times), "with didx_fit", didx_fits[np.argmax(engagement_times)%len(didx_fits)])
        if np.max(engagement_times) > -2e3:
            engagement_time = 0
        else:
            engagement_time = np.infty

    if plot:
        plt.plot(time_loadcell, load_cell, marker='o')

        for engagement_time_temp in engagement_times_filtered:
            didx_fit = didx_fits[np.where(engagement_times == engagement_time_temp)[0][0]]
            time_fit = time_loadcell[idx_relay_input_rise + didx_fit_start:idx_relay_input_rise + didx_fit]
            slope, intercept, r, p, se = linregress(time_fit, load_cell[idx_relay_input_rise + didx_fit_start:idx_relay_input_rise + didx_fit])
            time_fit = np.linspace((offset_when_voltage_enabled - intercept)/slope, time_fit[-1])
            plt.plot(time_fit, slope*time_fit + intercept, 'k', alpha=0.2)
        plt.plot(time_loadcell, relay_input, 'tab:red')
        plt.axhline(offset_when_voltage_enabled, color='k')
        # plt.plot(didx_fits, r2scores)
        plt.grid(True)
        plt.xlabel("Time (us)")
        plt.show()

    return engagement_time


def calculate_release_time(file_name, plot=False, verbose=True):
    time_loadcell, load_cell, relay_input, relay_output, idx_stop_moving, idx_relay_input_rise, idx_relay_input_fall = import_data(file_name,
                                                                                                                                   zero_load_cell=False,
                                                                                                                                   convert_to_kPa=False,
                                                                                                                                   convert_time_to_seconds=False)
    # The filtered version was just used to calculate the final converged value, and to check whether the load cell released unnaturally before the voltage released
    sos = butter(2, 250, btype='lowpass', fs=12800, output='sos')
    load_cell_filtered = sosfiltfilt(sos, load_cell)
    load_cell_final = np.nanmean(load_cell[idx_relay_input_fall + int(12800*0.01):idx_relay_input_fall + int(12800*0.11)])
    if np.where(load_cell_filtered > load_cell_final + 0.3*(np.max(load_cell_filtered) - load_cell_final))[0][-1] < idx_relay_input_fall:
        if verbose:
            print("Load cell released unnaturally before voltage released, released at",
                  time_loadcell[np.where(load_cell_filtered > load_cell_final + 0.3*(np.max(load_cell_filtered) - load_cell_final))[0][-1]])
        return np.inf, np.inf, np.inf, np.inf, np.inf

    load_cell_before_drop = np.max(load_cell[idx_relay_input_fall - int(12800*0.1):idx_relay_input_fall + int(12800*0.1)])
    loadcell_fall_threshold = load_cell_final + (load_cell_before_drop - load_cell_final)*0.1
    time_fall = time_loadcell[np.where((time_loadcell > time_loadcell[idx_relay_input_fall]) & (load_cell < loadcell_fall_threshold))[0][0]]
    release_time_ms = (time_fall - time_loadcell[idx_relay_input_fall])/1e3

    loadcell_fall_threshold_10 = load_cell_final + (load_cell_before_drop - load_cell_final)*0.9
    time_fall_10 = time_loadcell[np.where((time_loadcell > time_loadcell[idx_relay_input_fall]) & (load_cell < loadcell_fall_threshold_10))[0][0]]
    release_time_ms_10 = (time_fall_10 - time_loadcell[idx_relay_input_fall])/1e3

    if plot:
        plt.plot(time_loadcell, load_cell, marker='o')
        plt.plot(time_loadcell, load_cell_filtered, marker='o')
        plt.axvline(time_loadcell[idx_relay_input_fall], color='tab:red')
        plt.axhline(loadcell_fall_threshold, color='k')
        plt.axhline(load_cell_before_drop, color='k')
        plt.axhline(loadcell_fall_threshold_10, color='k')
        plt.axhline(load_cell_final, color='k')
        plt.grid(True)
        plt.xlabel("didx_fit")
        plt.show()

    load_cell_initial = np.nanmean(load_cell[idx_relay_input_rise - int(12800*0.005):idx_relay_input_rise])
    load_cell_max = np.nanmax(load_cell)
    return release_time_ms, load_cell_before_drop, load_cell_initial, load_cell_max, release_time_ms_10


if __name__ == '__main__':
    file_loc = "../data/"
    file_name = "20241027_17_12_31_bts_dx=1.5_v=0.5_Fpre=-1g_V=250_f=1000_w=2_brass"

    F_preload = calculate_Fpreload(file_loc + file_name + ".npy", w_pin=2e-3, depth_pin=2e-3)
    plot_engage_time = True
    engagement_time = calculate_engagement_time(file_loc + file_name + ".npy", plot=plot_engage_time)  # , verbose=True)
    release_time, load_cell_before_drop, load_cell_initial, load_cell_max, release_time_ms_10 = calculate_release_time(file_loc + file_name + ".npy", plot=not plot_engage_time)
    print("Engagement Time:", engagement_time, "us, Release Time:", release_time, "ms, Release Time 10%:", release_time_ms_10, ", F_preload:", F_preload, "N")
