#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os
import sys


def read_pairs(f):
    text_pairs = np.loadtxt(f, dtype='str', delimiter='\n')
    list_pairs = [[pair.split('|') for pair in line.strip().split(' ')] for line in text_pairs]
    return list_pairs


def distance(x0, x1, dimensions):
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))


def convert_pairs_to_tuples(pairs):
    tuples = []
    for pair in pairs:
        tuples.append((int(pair[0]), int(pair[1])))
    return tuples


def mode(arr):
    vals, counts = np.unique(arr, return_counts=True)
    index = np.argmax(counts)
    return vals[index], counts[index]


def get_polymer_connections(polymer_path, window=10, threshold=0.7):
    connections = []
    for i in range(len(polymer_path) - 2 * window):
        segment1 = np.array(polymer_path[i:i+window])
        segment2 = np.array(polymer_path[i+window:i+2*window])
        mode1, count1 = mode(segment1)
        mode2, count2 = mode(segment2)
        if mode1 != mode2 and count1 > window * threshold and count2 > window * threshold:
            connections.append((mode1, mode2))
            connections.append((mode2, mode1))
    return list(set(connections))


def get_capsomer_connected_state(capsomer, polymer_connections, original_bonds):
    num_connections = np.sum(np.array(polymer_connections)[:,0] == capsomer)
    if num_connections == 0:
        return 1
    if num_connections == 1:
        return 2
    if num_connections == 4:
        return 7
    if num_connections == 5:
        return 8
    connected_neighbors = np.array(polymer_connections)[np.array(np.array(polymer_connections)[:,0]==capsomer).T][:, 1]
    if num_connections == 2:
        if ((int(connected_neighbors[0]), int(connected_neighbors[1])) in original_bonds or
            (int(connected_neighbors[1]), int(connected_neighbors[0])) in original_bonds):
            return 3
        else:
            return 4
    for i in range(3):
        if ((int(connected_neighbors[i]), int(connected_neighbors[i-1])) not in original_bonds and
            (int(connected_neighbors[i-1]), int(connected_neighbors[i])) not in original_bonds and
            (int(connected_neighbors[i]), int(connected_neighbors[i-2])) not in original_bonds and
            (int(connected_neighbors[i-2]), int(connected_neighbors[i])) not in original_bonds):
            return 6
    return 5


def flatten(l):
    try:
        return flatten(l[0]) + (flatten(l[1:]) if len(l) > 1 else []) if type(l) is list else [l]
    except IndexError:
        return []


def write_pairs_to_file(filepath, pair_data):
    with open(filepath, 'w') as f:
        for frame in pair_data:
            line = ''
            for pair in frame:
                line = line + str(pair[0]) + '|' + str(pair[1]) + ' '
            f.write(line + '\n')
    return


deb = sys.argv[1]
att = sys.argv[2]
pss = sys.argv[3]
rep = sys.argv[4]

animation_period = 2
run_length = 5000
initial = 0
x_scale = np.arange(initial * animation_period, run_length, animation_period)
np.savetxt('./xscale/xscale.txt', x_scale)

pairs = read_pairs('./raw/pairs_DEB-%s_ATT-%s_PSS-%s_%s.txt' % (str(deb), str(att), str(pss), str(rep)))
bonds = np.loadtxt('./raw/bonds_DEB-%s_ATT-%s_PSS-%s_%s.txt' % (str(deb), str(att), str(pss), str(rep)))
polymer_paths = np.loadtxt('./raw/polymerpath_DEB-%s_ATT-%s_PSS-%s_%s.txt' % (str(deb), str(att), str(pss), str(rep)))

if not os.path.exists('./raw/polymerconnections2_DEB-%s_ATT-%s_PSS-%s_%s.txt' %
                      (str(deb), str(att), str(pss), str(rep))):
    polymer_trajectory = []
    for i in range(len(polymer_paths)):
        polymer_trajectory.append(get_polymer_connections(polymer_paths[i], window=7))
    write_pairs_to_file('./raw/polymerconnections2_DEB-%s_ATT-%s_PSS-%s_%s.txt' %
                        (str(deb), str(att), str(pss), str(rep)), polymer_trajectory)
else:
    polymer_trajectory = read_pairs('./raw/polymerconnections2_DEB-%s_ATT-%s_PSS-%s_%s.txt' %
                                    (str(deb), str(att), str(pss), str(rep)))

if not os.path.exists('./raw/capsomerstates2_DEB-%s_ATT-%s_PSS-%s_%s.txt' %
                      (str(deb), str(att), str(pss), str(rep))):
    capsomer_states = []
    # Identify original bonds
    original_pairs = convert_pairs_to_tuples(pairs[0])
    capsomers = list(set(np.array(original_pairs).flatten()))
    disassembly_index = len(bonds) - list(reversed(bonds)).index(30.) -  1
    for j in range(min(disassembly_index + 10, len(bonds))):
        capsomer_states.append([])
        # Identify polymer connections
        polymer_connections = polymer_trajectory[j]
        # Identify connected state for each capsomer
        for capsomer in capsomers:
            capsomer_states[-1].append(get_capsomer_connected_state(capsomer, polymer_connections, original_pairs))
    np.savetxt('./raw/capsomerstates2_DEB-%s_ATT-%s_PSS-%s_%s.txt' % (str(deb), str(att), str(pss), str(rep)),
               np.array(capsomer_states, dtype=int))

