#!/usr/bin/env python
# coding: utf-8

import gsd.hoomd
import numpy as np
import pandas as pd
import os
import sys

box_dim = np.array([40.0, 40.0, 40.0])


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


def get_polymer_path(traj, frame):
    particles = get_particles(traj, frame)
    p_particles = particles[particles['type']=='P']
    m_particles = particles[particles['type']=='M']
    nearest_bodies = []
    for i, p in p_particles.iterrows():
        nearest_body = m_particles['body'][575]
        nearest_distance = distance(np.array([p['position_x'], p['position_y'], p['position_z']]), 
                                    np.array([m_particles['position_x'][575], 
                                     m_particles['position_y'][575], 
                                     m_particles['position_z'][575]]), box_dim)
        for j, m in m_particles.iterrows():
            pm_distance = distance(np.array([p['position_x'], p['position_y'], p['position_z']]),
                                   np.array([m['position_x'], m['position_y'], m['position_z']]),
                                   box_dim)
            if pm_distance < nearest_distance:
                nearest_distance = pm_distance
                nearest_body = m['body']
        nearest_bodies.append(nearest_body)
    return nearest_bodies


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
p_paths = []

if not os.path.exists('./raw/polymerpath_DEB-%s_ATT-%s_PSS-%s_%s.txt' %
                      (str(deb), str(att), str(pss), str(rep))):
    for frame in range(initial, int(run_length / animation_period)):
        p_path = get_polymer_path(f, frame)
        p_paths.append(p_path)
    np.savetxt('./raw/polymerpath_DEB-%s_ATT-%s_PSS-%s_%s.txt' %
               (str(deb), str(att), str(pss), str(rep)), np.array(p_paths))

