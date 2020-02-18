import argparse
import numpy as np

def main(args):
    iname = args.iname
    coords = []
    for line in open(iname, 'rt'):
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        coords.append([x, y, z])
    coords = np.array(coords)
    center = [round(x, 3) for x in coords.mean(axis=0)]
    print(f'--center_x {x} --center_y {y} --center_z {z} --size_x 25 --size_y 25 --size_z 25')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('iname', type=str)
    args = parser.parse_args()
    main(args)
