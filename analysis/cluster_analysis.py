#!/usr/bin/env python
# coding: utf-8

import gsd.hoomd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys

box_dim = np.array([40.0, 40.0, 40.0])


def add_group(cluster, x, subunit_bonded_couples):
    new_bonds = 0
    for index, pair in enumerate(subunit_bonded_couples):
        if not pair:  # Couple has been removed previously
            continue
        elif pair[0] == cluster[x]:
            new_bonds += 1
            if pair[1] not in cluster:
                cluster.append(pair[1])
            subunit_bonded_couples[index] = []
        elif pair[1] == cluster[x]:
            new_bonds += 1
            if pair[0] not in cluster:
                cluster.append(pair[0])
            subunit_bonded_couples[index] = []
        else:
            continue
    return new_bonds


def find_clusters(sorted_attractor_coords, box):

    """
    Determines intercapsomer bonds and clusters of bonded capsomers

    :param sorted_attractor_coords: list of np.array, 'A' particle coordinates separated by capsomer/body
    :param box: list, 3D simulation box dimensions
    :return cluster_summary: list, cluster summary information - [[number of elements, number of bonds]...]
    """
    
    attachment_cutoff = .5
    attached_subunits = len(sorted_attractor_coords)
    subunit_bonded_couples = []  # All bonded pairs of subunits
    for i in range(attached_subunits):  # i, j are a distinct pair of subunits
        for j in range(0, i):
            attachment_counter = 0
            for i_A in sorted_attractor_coords[i]:  # i_A, j_A are attractors in i and j
                for j_A in sorted_attractor_coords[j]:
                    if not (j == i):
                        if distance(i_A, j_A, box) < attachment_cutoff:
                            attachment_counter += 1  # Marks one attachment
            if attachment_counter == 2:  # Two attachments means bond between subunits
                subunit_bonded_couples.append([i, j])

    clusters = []  # first element is number of bonds, second element is list of bonded subunits

    n = 0
    for pair in subunit_bonded_couples:
        iteration = 0
        if pair:  # Might have been removed by add_group in previous loops
            clusters.append([0, [pair[0]]])
            cluster_size = len(clusters[n][1])
            while iteration < cluster_size:
                # add_group adds subunits to the list of elements in the cluster
                new_bonds = add_group(clusters[n][1], iteration, subunit_bonded_couples)
                iteration += 1
                cluster_size = len(clusters[n][1])
                clusters[n][0] += new_bonds
            n += 1
    cluster_summary = []
    if not clusters:  # No clusters found
        clusters = [[0, []]]
    for cluster in clusters:
        cluster_summary.append([len(cluster[1]), cluster[0]])  # (number of subunits, number of bonds)
    cluster_summary.sort(reverse=True)  # Sort cluster list from largest to smallest
    return cluster_summary


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


def get_AA_contacts(traj, frame, center=None, radius=None):
    bond_length = 0.5
    particle_info = get_particles(traj, frame)
    A_particles = particle_info.loc[particle_info['type'] == 'A']
    A_coords = np.array([np.array(A_particles['position_x']),
                         np.array(A_particles['position_y']),
                         np.array(A_particles['position_z'])]).T
    if radius:
        A_coords = A_coords[np.where(distance(A_coords, center, box_dim) < radius)]
    distances = np.array([])
    for i in range(len(A_coords)):
        distances = np.append(distances, distance(A_coords, A_coords[i], box_dim))
    return (np.sum(distances < bond_length) - len(A_coords)) / 2


def get_PM_contacts(traj, frame, center=None, radius=None):  # Based on P-M proximity
    bond_length = 1
    particle_info = get_particles(traj, frame)
    P_particles = particle_info.loc[particle_info['type'] == 'P']
    M_particles = particle_info.loc[particle_info['type'] == 'M']
    P_coords = np.array([P_particles['position_x'],
                         P_particles['position_y'],
                         P_particles['position_z']]).T
    M_coords = np.array([M_particles['position_x'],
                         M_particles['position_y'],
                         M_particles['position_z']]).T
    if radius:  # For efficiency
        M_coords = M_coords[np.where(distance(M_coords, center, box_dim) < radius)]
    adsorbed = 0
    for i in range(len(M_coords)):
        distances = distance(P_coords, M_coords[i], box_dim)
        if np.any(distances < bond_length):
            adsorbed += 1
    return adsorbed


def get_subs_in_capsid(traj, frame, center=None, radius=None):
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

    return find_clusters(A_coord_select, box_dim)


deb = sys.argv[1]
att = sys.argv[2]
pss = sys.argv[3]
rep = sys.argv[4]

animation_period = 2
run_length = 5000
initial = 0
x_scale = np.arange(initial * animation_period, run_length, animation_period)
np.savetxt('./xscale/xscale.txt', x_scale)

adsorbed_subunits = []
max_cluster = []
cluster_bonds = []

f = gsd.hoomd.open('../DEB-%s_ATT-%s_PSS-%s_%s/traj.gsd' % (str(deb), str(att), str(pss), str(rep)))

if not os.path.exists('./raw/bonds_DEB-%s_ATT-%s_PSS-%s_%s.txt' % (str(deb), str(att), str(pss), str(rep))):
    for frame in range(initial, int(run_length / animation_period)):
        P_center = get_P_center(f, frame)
        adsorbed_subunits.append(get_PM_contacts(f, frame, center=P_center, radius=40))
        cluster_data = get_subs_in_capsid(f, frame, center=P_center, radius=40)
        max_cluster.append(cluster_data[0][0])
        cluster_bonds.append(cluster_data[0][1])
    np.savetxt('./raw/adsorbed_DEB-%s_ATT-%s_PSS-%s_%s.txt' % (str(deb), str(att), str(pss), str(rep)), adsorbed_subunits)
    np.savetxt('./raw/cluster_DEB-%s_ATT-%s_PSS-%s_%s.txt' % (str(deb), str(att), str(pss), str(rep)), max_cluster)
    np.savetxt('./raw/bonds_DEB-%s_ATT-%s_PSS-%s_%s.txt' % (str(deb), str(att), str(pss), str(rep)), cluster_bonds)

