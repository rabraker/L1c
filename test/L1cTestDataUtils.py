#!/usr/bin/env python3
import numpy as np
import json
import codecs


def jsonify(data):

    for key, value in data.items():
        if type(value) is np.ndarray:
            if len(value.flatten()) == 1:
                data[key] = value.flatten()[0]
            else:
                data[key] = value.flatten().tolist()

    return data


def save_json(data, file_path):

    with codecs.open(file_path, 'w') as fid:
        json.dump(data, fid,  separators=(',', ':'), sort_keys=True, indent=4)

    # json.dump(data, codecs.open(file_path, 'w'),
    #           separators=(',', ':'), sort_keys=True, indent=4)


def data_dir():
    import os

    data_dir_parent = os.getenv("srcdir")
    if data_dir_parent is None:
        data_dir_parent = "."

    return data_dir_parent+"/test_data"
