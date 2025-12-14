from parse_data_filenames import *
from calculate_timing_stats import *
from datetime import datetime

save_folder = "../data/"
# file_name = "20241206_23_17_24_timing_alldata"
file_name = None

if file_name is None:
    now = datetime.now()
    name_clarifier = "_timing_alldata"
    timestamp = now.strftime("%Y%m%d_%H_%M_%S") + name_clarifier
    save_folder = "../data/"

    file_loc = "../data/"
    all_file_names = get_file_dict("20241027_15_25_21", "20241107_15_29_44", key='V', verbose=False)

    all_V = set()
    for V, file_names in all_file_names.items():
        for file_name, speed, V, freq, Fpreload, width in file_names:
            all_V.add(V)
    all_V = sorted(all_V)

    num_files = 0
    num_skipped = 0
    all_engagement_times = []
    all_release_times = []
    all_release_times_10pct = []
    all_data = []
    for V, file_names in all_file_names.items():
        print("V =", V, len(all_file_names[V]))
        for file_name, speed, V, freq, Fpreload, width in file_names:
            print('----------', file_name, '----------')
            depth_pin = get_pin_depth_from_width(width)
            engagement_time = calculate_engagement_time(file_loc + file_name)
            release_time, load_cell_before_drop, load_cell_initial, load_cell_max, release_time_10pct = calculate_release_time(file_loc + file_name)
            if engagement_time > 50:
                print("Skipping, engagement time unnaturally large:", engagement_time, 'Release Time', release_time)
                num_skipped += 1
                continue
            if release_time > 10:
                print("Skipping, release time unnaturally large:", release_time, 'Engagement Time', engagement_time)
                num_skipped += 1
                continue
            num_files += 1
            all_engagement_times.append(engagement_time)
            all_release_times.append(release_time)
            all_release_times_10pct.append(release_time_10pct)
            all_data.append((file_name, speed, V, freq, Fpreload, width, engagement_time, release_time, load_cell_before_drop, load_cell_initial, load_cell_max, release_time_10pct))
    print("Total # of Files:", num_files, "Num Skipped", num_skipped)
    print("Stats on Engagement Times:")
    print("Mean", np.nanmean(all_engagement_times), "Std", np.nanstd(all_engagement_times), "Median", np.nanmedian(all_engagement_times))
    print(all_engagement_times)

    print("Stats on Release Times:")
    print("Mean", np.nanmean(all_release_times), "Std", np.nanstd(all_release_times), "Median", np.nanmedian(all_release_times))
    print(all_release_times)

    print("Stats on Release Times 10pct:")
    print("Mean", np.nanmean(all_release_times_10pct), "Std", np.nanstd(all_release_times_10pct), "Median", np.nanmedian(all_release_times_10pct))
    print(all_release_times_10pct)

    # np.save(save_folder + timestamp + "_tr10pct", [all_V, all_engagement_times, all_release_times, all_release_times_10pct, all_data])
else:
    all_V, all_engagement_times, all_release_times, all_release_times_10pct, all_data = np.load(save_folder + file_name + ".npy", allow_pickle=True)

    # print("# Total:", len())
    all_engagement_times = np.array(all_engagement_times)
    all_release_times = np.array(all_release_times)

    num_files = len(all_engagement_times)
    print("Total # of Files:", num_files)
    print("Stats on Engagement Times:")
    print("Mean", np.nanmean(all_engagement_times), "Std", np.nanstd(all_engagement_times), "Median", np.nanmedian(all_engagement_times))
    print("# Engagement Times < 25us", len(all_engagement_times[all_engagement_times <= 25]), len(all_engagement_times[all_engagement_times <= 25])/num_files, '%')
    print("# Engagement Times < 15us", len(all_engagement_times[all_engagement_times <= 15]), len(all_engagement_times[all_engagement_times <= 15])/num_files, '%')
    print(list(all_engagement_times))

    print("Stats on Release Times:")
    print("Mean", np.nanmean(all_release_times), "Std", np.nanstd(all_release_times), "Median", np.nanmedian(all_release_times))
    print("# Release Times < 1 ms", len(all_engagement_times[all_release_times <= 1]), len(all_release_times[all_release_times <= 1])/num_files, '%')
    print(list(all_release_times))

    file_name = "20241214_05_46_46_timing_alldata_tr10pct"
    _, _, _, all_release_times_10pct, _ = np.load(save_folder + file_name + ".npy", allow_pickle=True)
    all_release_times_10pct = np.array(all_release_times_10pct)
    print("Stats on Release Times 10pct:")
    print("Mean", np.nanmean(all_release_times_10pct), "Std", np.nanstd(all_release_times_10pct), "Median", np.nanmedian(all_release_times_10pct))
    print(list(all_release_times_10pct))
    print("Average delay between 90% and 10% release times", np.mean(all_release_times - all_release_times_10pct), np.std(all_release_times - all_release_times_10pct))
