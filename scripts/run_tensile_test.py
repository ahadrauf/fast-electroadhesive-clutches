"""
Run a tensile test. The code in the teensy/tensile_tester_motor_drive_loadcell_logdata folder should be loaded onto the Teensy.
"""
import numpy as np
import csv
import matplotlib
import matplotlib.pyplot as plt
import serial
import threading
from datetime import datetime, timedelta
import logging
import time
import struct
from scipy.stats import linregress
import beepy

# Teensy settings
teensy_port = "COM5"
teensy = serial.Serial(teensy_port)
teensy.flushInput()
time.sleep(1)
num_runs = 1

for idx_run in range(num_runs):
    now = datetime.now()
    name_clarifier = "_bts_"
    timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
    print(timestamp)

    distance_brake_engaged = 1.5  # 6  # 4  # mm
    speed_brake_engaged = 0.5  # mm/s
    V = 200
    f = 20000
    w_pin = 2
    depth_pin = 2
    clarifier = "dx={}_v={}_Fpre=-1g_V={}_f={}_w={}_brass".format(distance_brake_engaged, speed_brake_engaged, V, f, w_pin)

    distance_multiplier = 0.1  # copied from Arduino code
    speed_multiplier = 0.05
    NEMA14_STEPS_PER_REV = 200
    NEMA14_LEAD_SCREW_PITCH = 2
    NEMA14_MICROSTEPS = 256
    NEMA14_STEPS_PER_MM = NEMA14_STEPS_PER_REV * NEMA14_MICROSTEPS / NEMA14_LEAD_SCREW_PITCH
    DELAY_STEP = int(500000 / (int(speed_brake_engaged / speed_multiplier) * speed_multiplier) / NEMA14_STEPS_PER_MM)

    # print("Actual speed:", 500000 / DELAY_STEP / NEMA14_STEPS_PER_MM, "Delay step", DELAY_STEP)

    # Logging settings
    save_every = 1
    terminate_lock = threading.Lock()
    run_count = 1
    show_every_measurement = 500
    average_run_time = 1.0
    stop_running = False
    count_input_errors = 0


    def save_data(time_data, load_cell_data, moving_forward_data, relay_output_data, clarifier=""):
        csv_file_path = '../data/strain_tests/' + timestamp + clarifier + '.csv'
        all_data = np.stack([time_data, load_cell_data, moving_forward_data, relay_output_data]).transpose()
        with open(csv_file_path, "w", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Time", "Load Cell", "Moving Forward", "Relay Output"])
            writer.writerows(all_data)
        np.save('../data/strain_tests/' + timestamp + clarifier + '.npy', all_data, allow_pickle=True)


    def read_data():
        global count_input_errors
        ser_bytes = teensy.readline()
        decoded_str = ser_bytes[0:len(ser_bytes) - 2].decode("utf-8")
        if len(decoded_str) == 0:
            return None
        decoded_str_split_by_decimal = decoded_str.split('\r')
        if len(decoded_str_split_by_decimal) > 2:
            print("Error in decoding string:", decoded_str, flush=True)
            decoded_str = decoded_str_split_by_decimal[0]  # + '.' + decoded_str_split_by_decimal[1]
        try:
            decoded_bytes = float(decoded_str)
            return decoded_bytes
        except Exception as e:
            print("WARNING: ", e, decoded_str, decoded_str_split_by_decimal)
            logging.warning(e)
            count_input_errors += 1
            return None


    def plot_data(time_data, load_cell_data, moving_forward_data, relay_output_data):
        # Error correction on time data
        for i in range(100, len(time_data) - 100):
            if time_data[i] == np.max(time_data[i - 100: i + 100]) or time_data[i] == np.min(time_data[i - 100: i + 100]):
                time_data[i] = (time_data[i - 1] + time_data[i + 1]) / 2

        fig, ax1_right = plt.subplots(1, 1, layout='constrained')
        ax1 = fig.add_subplot(111, sharex=ax1_right, frameon=False)
        ax1.plot(time_data / 1e6, load_cell_data, c='tab:blue')
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Load Cell Reading (N)')

        ax1_right.plot(time_data / 1e6, relay_output_data, 'tab:red', alpha=0.75)
        ax1_right.plot(time_data / 1e6, moving_forward_data * 100, 'tab:green', alpha=0.75)
        ax1_right.yaxis.tick_right()
        ax1_right.yaxis.set_label_position("right")
        plt.ylabel("Voltage (V)")
        ax1_right.yaxis.label.set_color('tab:red')
        ax1_right.tick_params(axis='y', colors='tab:red')

        ax1.grid()
        plt_title = timestamp + clarifier
        if plt_title is not None:
            fig.suptitle(plt_title)
        plt.savefig("../figures/strain_tests/" + timestamp + clarifier + ".png")
        if num_runs == 1:
            plt.show()
        else:
            plt.close(fig)


    def print_stats(time_data, load_cell_data, moving_forward_data, relay_output_data):
        print("Average measurement time", np.mean(time_data[1:] - time_data[:-1]) * 1e-3, "ms, Frequency:",
              1e6 / np.mean(time_data[1:] - time_data[:-1]), "Hz")
        print("Initial force:", np.nanmean(load_cell_data[:9000]))
        print("Max force:", np.nanmax(load_cell_data))
        print("dF:", np.nanmax(load_cell_data) - np.nanmean(load_cell_data[:9000]))
        didx_fit = int(12800 * 0.003)
        idx_relay_input_rise = max(np.where(relay_output_data > 50)[0][0], np.where(moving_forward_data > 0.5)[0][0])
        loadcell_fit = load_cell_data - np.nanmean(load_cell_data[idx_relay_input_rise - didx_fit * 4:idx_relay_input_rise])
        idx_relay_output_fall = max(np.where(relay_output_data > 50)[0][-1] + 1, np.where(moving_forward_data > 0.5)[0][-1] + 1)
        slope, intercept, r, p, se = linregress(time_data[idx_relay_input_rise:idx_relay_input_rise + didx_fit],
                                                loadcell_fit[idx_relay_input_rise:idx_relay_input_rise + didx_fit])
        print("Voltage rose at t =", time_data[idx_relay_input_rise], 'Voltage fell at t=', time_data[idx_relay_output_fall])
        print("Engagement time", (-intercept / slope - time_data[idx_relay_input_rise]) / 1e3, "ms")

        mu_objet_s = 0.281  # 0.173  # 0.28
        mu_pvdf_s = 0.188  # 0.154  # 0.38
        mu_objet_k = 0.173
        mu_pvdf_k = 0.154
        mu_tot_k = mu_objet_k + mu_pvdf_k
        mu_tot_s = mu_objet_s + mu_pvdf_s
        initial_friction = np.nanmean(load_cell_data[idx_relay_input_rise - int(12800 * 0.1):idx_relay_input_rise])
        m_pin = 8.049e3 * w_pin * depth_pin * 1e-6 * 100e-3
        m_pvdf = 1.78e3 * 24e-6 * 60e-3 * w_pin
        g = 9.81
        F_preload = (1 / mu_tot_k) * (initial_friction - m_pin * g * mu_objet_k - m_pvdf * g * mu_tot_k)
        print("Preload force:", F_preload, "Initial friction:", initial_friction)

        load_cell_final = np.nanmean(load_cell_data[-1000:])
        load_cell_before_drop = load_cell_data[idx_relay_output_fall]
        loadcell_fall_threshold = load_cell_final + (load_cell_before_drop - load_cell_final) * 0.1
        time_fall = time_data[np.where((time_data > time_data[idx_relay_output_fall]) & (load_cell_data < loadcell_fall_threshold))[0][0]]
        print("Fall time", (time_fall - time_data[idx_relay_output_fall]) / 1e3, "ms", time_data[idx_relay_output_fall], time_fall, "threshold", loadcell_fall_threshold)
        print("Average sample frequency while voltage is on", 1e6 / np.mean(time_data[idx_relay_input_rise + 1:idx_relay_output_fall] -
                                                                            time_data[idx_relay_input_rise:idx_relay_output_fall - 1]), "Hz")


    base_clarifier = clarifier
    orig_start_time = time.perf_counter_ns()

    ##############################################################################################
    # Run Actuonix motor
    ##############################################################################################
    try:
        start_time = time.perf_counter_ns()
        teensy.write(struct.pack('>BB', int(distance_brake_engaged / distance_multiplier),
                                 int(speed_brake_engaged / speed_multiplier)))

        start_reading_time = time.perf_counter_ns()
        num_measurements = read_data()
        if num_measurements is not None:
            print("Num measurements:", num_measurements)
            num_measurements = int(num_measurements)
        else:
            print("ERROR: couldn't read the number of measurements")
        time_data = np.zeros(num_measurements)
        load_cell_data = np.zeros(num_measurements)
        moving_forward_data = np.zeros(num_measurements)
        relay_output_data = np.zeros(num_measurements)
        # beepy.beep(sound='coin')

        curr_data = read_data()
        print("Reading time data, curr_data =", curr_data, flush=True, end=" ")
        curr_idx = 0
        while curr_data != num_measurements:
            if curr_idx % 50000 == 0:
                print(curr_idx, end=" | ", flush=True)
            time_data[curr_idx] = curr_data
            curr_data = read_data()
            curr_idx += 1

        curr_data = read_data()
        print("\nReading load cell data, curr_data =", curr_data, "last time =", time_data[-1], ",", flush=True, end=" ")
        curr_idx = 0
        while curr_data != num_measurements:
            if curr_idx % 50000 == 0:
                print(curr_idx, end=" | ", flush=True)
            load_cell_data[curr_idx] = curr_data
            if curr_data is not None and abs(curr_data) < 1000:
                load_cell_data[curr_idx] = curr_data
            else:
                load_cell_data[curr_idx] = np.nan
            curr_data = read_data()
            curr_idx += 1

        curr_data = read_data()
        print("\nReading relay input data, curr_data =", curr_data, "max_force =", np.nanmax(load_cell_data), flush=True, end=" ")
        curr_idx = 0
        while curr_data != num_measurements:
            if curr_idx % 50000 == 0:
                print(curr_idx, end=" | ", flush=True)
            moving_forward_data[curr_idx] = curr_data
            curr_data = read_data()
            curr_idx += 1

        curr_data = read_data()
        print("\nReading relay output data, curr_data =", curr_data, flush=True, end=" ")
        curr_idx = 0
        while curr_data != num_measurements:
            if curr_idx % 50000 == 0:
                print(curr_idx, end=" | ", flush=True)
            relay_output_data[curr_idx] = curr_data
            curr_data = read_data()
            curr_idx += 1
        print("Max voltage read:", np.nanmax(relay_output_data))

        end_time = time.perf_counter_ns()

        stop_running = True
        save_data(time_data, load_cell_data, moving_forward_data, relay_output_data, clarifier=clarifier)
        plot_data(time_data, load_cell_data, moving_forward_data, relay_output_data)

        print('\n---------------------Run {}--------------------'.format(idx_run + 1))
        print("Total time to read data:", (end_time - start_reading_time) / 1e9, "Avg. time/reading:",
              (end_time - start_reading_time) / 1e6 / (4 * num_measurements), "ms")
        print('-----------------------------------------------')
    except KeyboardInterrupt as e:
        print(e)
        terminate_lock.acquire()
        terminate_lock.release()
        stop_running = True
    else:
        terminate_lock.acquire()
        terminate_lock.release()
        run_count += 1

    # time.sleep(3)
    # beepy.beep(sound='ready')
    time.sleep(2)
