import os
from pathlib import PurePath
from glob import glob
import shutil
import json
import argparse

import mdtraj
from sklearn.cluster import AgglomerativeClustering
import numpy as np
from matplotlib import pyplot as plt


def main1(base_path, prefix='rank', cut=50.0):
    print('Clustering Trajectories...')
    print(f'  {base_path}')

    prmtop = os.path.join(base_path, 'init.prmtop')
    traj = glob(os.path.join(base_path, 'prod.?.xtc'))
    traj.sort()

    print(traj)

    univ = mdtraj.load(traj, top=prmtop)
    mask = univ.topology.select('protein and not type "H"')
    prot = univ.atom_slice(mask)

    if not os.path.exists(f'{prefix}.dist.npy'):
        nn = len(prot)
        distance_matrix = np.zeros((nn, nn))

        print('  Calculate a distance matrix')
        for n in range(nn):
            distance_matrix[n] = mdtraj.rmsd(prot, prot, frame=n)
        np.save(f'{prefix}.dist.npy', distance_matrix)

    dist_mat = np.load(f'{prefix}.dist.npy')

    # plt.figure(figsize=(8, 8))
    # plt.imshow(dist_mat)
    # plt.show()

    print('  Clustering...')
    # AC dist_threshold --> 50.0
    clustering = AgglomerativeClustering(n_clusters=None, affinity='euclidean',
                                         compute_full_tree=True,
                                         distance_threshold=cut).fit(dist_mat)
    centers = []
    for label in set(clustering.labels_):
        l_mask = clustering.labels_ == label
        num_label = np.count_nonzero(l_mask)
        index_map = np.arange(0, len(clustering.labels_))[l_mask]
        block_mat = dist_mat[index_map, :][:, index_map]

        # print(index_map)
        # print(block_mat)

        # med = np.argmin(block_mat.sum(axis=0))  # L1 norm
        med = np.argmin((block_mat ** 2).sum(axis=0))  # L2 norm
        stat = {
            'average_distance': np.mean(block_mat[med]),
            'std': np.sqrt(np.mean(block_mat[med] ** 2))
        }
        # print(label, '-->', num_label, '[', index_map[med], ']')
        centers.append((label, num_label, index_map[med], stat))

    print('  Save...')
    centers.sort(key=lambda x: x[1], reverse=True)
    ranking_info = []
    for rank, (_, num, index, stat) in enumerate(centers):
        frame = prot[index]
        # frame.save_pdb(f'{prefix}.pdb')
        out = f'{prefix}_{rank:02d}.pdb'
        frame.save_pdb(out)
        stat = {k: float(v) for k, v in stat.items()}
        ranking_info.append(
            {'rank': int(rank), 'population': int(num), 'frame': int(index), 'stat': stat, 'output': str(out)})
        print(rank, '-->', num, '[', index, ']')
    with open(f'{prefix}.json', 'w') as f:
        json.dump({
            'prefix': prefix,
            'distance_threshold': cut,
            'ranking_info': ranking_info
        }, f, indent=4, sort_keys=False)


def main2():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--target-dir')
    parser.add_argument('-cut', '--threshold', default=50.0, type=float)
    args = parser.parse_args()
    targets = glob(os.path.join(args.target_dir, '*/'))
#    ranking = {'model_1': 2, 'model_2': 4, 'model_3': 1, 'model_4': 5, 'model_5': 3}  # H1106
    ranking = {}
    for target in targets:
        pp = PurePath(target)
        name = pp.parts[-1]
        print('=== Clustering ===')
        print(pp, '-->', name)
        try:
           main1(os.path.abspath(pp), name, args.threshold)
        except Exception as e:
            print(e)
    for name, rank in ranking.items():
        shutil.copy(f'{name}_00.pdb', f'r{rank}.pdb')
    with open('rank.json', 'w') as fp:
        json.dump(ranking, fp, indent=4)


if __name__ == '__main__':
    main2()
