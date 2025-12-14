import numpy as np


def import_data(file_loc, width_pattern=2, length_pattern=55.5, zero_load_cell=True, convert_to_kPa=True, convert_time_to_seconds=True):
    data = np.load(file_loc)
    has_relay_output = (data.shape[1] == 4)
    time_loadcell, load_cell, stage = data[:, 0], data[:, 1], data[:, 2]
    if has_relay_output:
        relay_output = data[:, 3]
    else:
        relay_output = None

    if convert_to_kPa:
        load_cell = load_cell/width_pattern/length_pattern*1000  # convert from N to kPa

    # Delete NaN values
    to_delete = []
    for idx in range(len(time_loadcell)):
        if np.isnan(time_loadcell[idx]) or np.isnan(load_cell[idx]) or load_cell[idx] > 50 or time_loadcell[idx] > 1e8:
            to_delete.append(idx)
        if idx >= 1 and abs(time_loadcell[idx] - time_loadcell[idx - 1]) >= 1000:
            # to_delete.append(idx)
            if time_loadcell[-1] == 0:  # there was a delay reading the time values, so some values got skipped
                idx_start_zeros = np.where(time_loadcell == 0)[0][0]
                num_nonzero_readings = idx_start_zeros - idx
                time_loadcell[-num_nonzero_readings:] = time_loadcell[idx:idx_start_zeros]
                time_loadcell[idx:len(time_loadcell) - num_nonzero_readings] = np.linspace(time_loadcell[idx - 1], time_loadcell[-num_nonzero_readings], len(time_loadcell) - num_nonzero_readings - (idx - 1) + 1)[1:-1]
                print("Moving back {} readings to zero range, starting at idx={}".format(num_nonzero_readings, idx_start_zeros))
            else:
                to_delete.append(idx)
        elif idx >= 1 and time_loadcell[idx] <= time_loadcell[idx - 1]:
            to_delete.append(idx)
    print("Deleted {} points".format(len(to_delete)))
    time_loadcell = np.delete(time_loadcell, to_delete)
    load_cell = np.delete(load_cell, to_delete)
    stage = np.delete(stage, to_delete)
    if has_relay_output:
        relay_output = np.delete(relay_output, to_delete)
    if convert_time_to_seconds:
        time_loadcell /= 1e6

    if has_relay_output:
        if np.max(relay_output) > 50:
            idx_relay_input_rise = max(np.where(stage > 0.5)[0][0], np.where(relay_output > 50)[0][0])
            idx_relay_input_fall = max(np.where(stage > 0.5)[0][-1], np.where(relay_output > 50)[0][-1]) + 1
        else:
            idx_relay_input_rise = np.where(stage > 0.5)[0][0]
            idx_relay_input_fall = np.where(stage > 0.5)[0][-1] + 1
        idx_stop_moving = idx_relay_input_fall
    else:
        idx_relay_input_rise = 0
        idx_stop_moving = 0
        idx_relay_input_fall = 0
        for i in range(1, len(time_loadcell)):
            if stage[i - 1] < 0.5 < stage[i]:
                if idx_relay_input_rise == 0:
                    idx_relay_input_rise = i
            if stage[i] < 0.5 < stage[i - 1]:
                if idx_stop_moving == 0:
                    idx_stop_moving = i
                else:
                    idx_relay_input_fall = i
    if zero_load_cell:
        load_cell -= np.nanmean(load_cell[idx_relay_input_rise - 100:idx_relay_input_rise - 1])

    return time_loadcell, load_cell, stage, relay_output, idx_stop_moving, idx_relay_input_rise, idx_relay_input_fall
