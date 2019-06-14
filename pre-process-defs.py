import re

# import ipdb
import sys

defs = {'NESTA_VERBOSE_IDX': 0,
        'NESTA_SIGMA_IDX': 1,
        'NESTA_MU_IDX': 2,
        'NESTA_TOL_IDX': 3,
        'NESTA_ALPV_IDX': 4,
        'NESTA_ALPH_IDX': 5,
        'NESTA_NCONT_IDX': 6,
        'NESTA_MODE_IDX': 7,
        'NESTA_DCTMODE_IDX': 8}


def do_sub(defs, line_, token):
    # keys =
    for key in defs.keys():
        # print(line)
        repl = '%d' % defs[key]
        line_ = re.sub(token+key+token, repl, line_)

    return line_


def compile_file(fpath_in, fpath_out):
    fid_in = open(fpath_in, 'r')
    fid_out = open(fpath_out, 'w+')
    token = '@'
    while True:
        line = fid_in.readline()
        if line == '':
            break

        line = do_sub(defs, line, token)
        fid_out.writelines(line)

    fid_in.close()
    fid_out.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise Exception("%s takes exactly two arguments"%sys.argv[0])

    fpathout = sys.argv[1]
    fpathin = sys.argv[2]

    compile_file(fpathin, fpathout)
