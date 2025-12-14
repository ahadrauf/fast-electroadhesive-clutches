from os import listdir
from os.path import isfile, join
import re


def get_file_dict(min_date, max_date, key, Fpreload_filter=None, V_filter=None, f_filter=None, width_filter=None, verbose=False):
    # file_loc = "../data/strain_tests/"
    file_loc = "../data/"
    float_template = r"([+-]?([0-9]+([.][0-9]*)?|[.][0-9]+))"
    filename_template = r"(\d+_\d\d_\d\d_\d\d)_bts_dx=" + float_template + r"_v=" + float_template + \
                        r"_Fpre=" + float_template + r"g_V=(\d*)_f=" + float_template + "_w=" + \
                        float_template + "_brass"
    banned_file_list = ["20241103_17_16_56_bts_dx=1.5_v=0.5_Fpre=-1g_V=150_f=1000_w=2_brass",
                        "20241029_15_31_38_bts_dx=1.5_v=0.5_Fpre=-1g_V=250_f=0.1_w=2_brass"]
    banned_file_list = [bf + '.npy' for bf in banned_file_list]

    all_data_files = [f for f in listdir(file_loc) if isfile(join(file_loc, f)) and f[-4:] == '.npy']
    data_files = {}
    for data_file in all_data_files:
        match = re.search(filename_template, data_file)
        if match and min_date <= match.group(1) <= max_date:
            date, dx, speed, Fpreload = match.group(1), float(match.group(2)), float(match.group(5)), float(match.group(8))
            V, f, width = float(match.group(11)), float(match.group(12)), float(match.group(15))
            if Fpreload_filter and Fpreload != Fpreload_filter:
                continue
            if V_filter and V != V_filter:
                continue
            if f_filter and f != f_filter:
                continue
            if width_filter and width != width_filter:
                continue
            if data_file in banned_file_list:
                continue
            data_files[data_file] = {'date': date, 'dx': dx, 'speed': speed, 'Fpreload': Fpreload, 'V': V, 'f': f, 'width': width}

    # Reorganize data for the desired variable
    # var = 'f'
    output = {}
    for filename in data_files.keys():
        data_file = data_files[filename]
        if data_file[key] not in output:
            output[data_file[key]] = []
        output[data_file[key]] += [(filename, data_file['speed'], data_file['V'], data_file['f'], data_file['Fpreload'], data_file['width'])]

    if verbose:
        print(output)
        print('{', end='')
        # Format files for insertion into another Python script
        for key in sorted(output.keys()):
            print(key, end=': [')
            for data_file in output[key]:
                if data_file == output[key][-1]:
                    print(str(data_file), '],')
                else:
                    print(str(data_file), ',')
        print('}')
    return output


def get_pin_depth_from_width(pin_width):
    if pin_width == 2 or pin_width == 4:
        return 2
    elif pin_width == 2.5 or pin_width == 5:
        return 2.5
    elif pin_width == 3 or pin_width == 6:
        return 3
    else:
        return 0


if __name__ == '__main__':
    get_file_dict("20241023_13_27_26", "20241023_14_54_16", key='Fpreload', verbose=True)
