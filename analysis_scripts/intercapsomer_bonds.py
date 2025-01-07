#!/usr/bin/env python
# coding: utf-8

import gsd.hoomd
import numpy as np
import pandas as pd
import os
import sys

box_dim = np.array([40.0, 40.0, 40.0])


def find_bonded_pairs(sorted_attractor_coords, rigid_bodies, box):

    """
    Determines intercapsomer bonded pairs

    :param sorted_attractor_coords: list of np.array, 'A' particle coordinates separated by capsomer/body
    :param rigid_bodies: list, rigid body numbers corresponding to the sorted list
    :param box: list, 3D simulation box dimensions
    :return subunit_bonded_pairs: list, bonded pairs of capsomers identified by body number
    """
    
    attachment_cutoff = .5
    attached_subunits = len(sorted_attractor_coords)
    bonded_pairs = []  # All bonded pairs of subunits
    for i in range(attached_subunits):  # i, j are a distinct pair of subunits
        for j in range(0, i):
            attachment_counter = 0
            for i_A in sorted_attractor_coords[i]:  # i_A, j_A are attractors in i and j
                for j_A in sorted_attractor_coords[j]:
                    if not (j == i):
                        if distance(i_A, j_A, box) < attachment_cutoff:
                            attachment_counter += 1  # Marks one attachment
            if attachment_counter == 2:  # Two attachments means bond between subunits
                bonded_pairs.append([rigid_bodies[i], rigid_bodies[j]])
    return bonded_pairs


def get_particles(traj, frame):
    particle_info = {
        'type': [f.read_frame(frame).particles.types[typeid]
                 for typeid in traj.read_frame(frame).particles.typeid],
        'body': f.read_frame(frame).particles.body,
        'position_x': traj.read_frame(frame).particles.position[:, 0],
        'position_y': traj.read_frame(frame).particles.position[:, 1],
        'position_z': traj.read_frame(frame).particles.position[:, 2],
    }
    return pd.DataFrame(particle_info)


def distance(x0, x1, dimensions):
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))


def get_P_center(traj, frame):
    particle_info = get_particles(traj, frame)
    P_particles = particle_info.loc[particle_info['type'] == 'P']
    if len(P_particles) == 0:
        return [[0.0, 0.0, 0.0]]
    P_coords = np.array([P_particles['position_x'],
                         P_particles['position_y'],
                         P_particles['position_z']]).T
    P_shift = P_coords[0]
    P_shifted = ((P_coords - P_shift) + (box_dim / 2)) % box_dim
    P_mean = (np.mean(P_shifted, axis=0).reshape(1, 3) + P_shift) % box_dim - (box_dim / 2)
    return P_mean


def get_pairs_in_capsid(traj, frame, center=None, radius=None):
    particle_info = get_particles(traj, frame)
    A_particles = particle_info.loc[particle_info['type'] == 'A']
    A_coords = np.array([np.array(A_particles['position_x']),
                         np.array(A_particles['position_y']),
                         np.array(A_particles['position_z'])]).T
    # Select particles within radius of center
    if radius:
        A_select = A_particles.iloc[np.where(distance(A_coords, center, box_dim) < radius)]
    else:
        A_select = A_particles
        
    A_bodies = A_select['body'].unique()
    A_coord_select = []
    for body in A_bodies:
        A_body_select = A_particles.loc[A_particles['body'] == body]  # A particles belonging to the same rigid body
        assert (len(A_body_select['body']) == 5)                      # There should be 5 of them
        A_coord_select.append(np.array([np.array(A_body_select['position_x']),
                                        np.array(A_body_select['position_y']),
                                        np.array(A_body_select['position_z'])]).T)

    return find_bonded_pairs(A_coord_select, A_bodies, box_dim)


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

f = gsd.hoomd.open('../DEB-%s_ATT-%s_PSS-%s_%s/traj.gsd' % (str(deb), str(att), str(pss), str(rep)))
pairs = []

if not os.path.exists('./raw/pairs_DEB-%s_ATT-%s_PSS-%s_%s.txt' %
                      (str(deb), str(att), str(pss), str(rep))):
    for frame in range(initial, int(run_length / animation_period)):
        P_center = get_P_center(f, frame)
        frame_bonded_pairs = get_pairs_in_capsid(f, frame, center=P_center, radius=40)
        pairs.append(frame_bonded_pairs)
    write_pairs_to_file('./raw/pairs_DEB-%s_ATT-%s_PSS-%s_%s.txt' %
                       (str(deb), str(att), str(pss), str(rep)), pairs)
