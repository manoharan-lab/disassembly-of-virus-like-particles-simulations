#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from config_dict import ConfigDict
import gsd.hoomd
import sys


def get_particles(traj, frame):
    particle_info = {
        'mass': traj.read_frame(frame).particles.mass,
        'diameter': traj.read_frame(frame).particles.diameter,
        'type': [traj.read_frame(frame).particles.types[typeid] 
                 for typeid in traj.read_frame(frame).particles.typeid],
        'body': traj.read_frame(frame).particles.body,
        'position_x': traj.read_frame(frame).particles.position[:, 0],
        'position_y': traj.read_frame(frame).particles.position[:, 1],
        'position_z': traj.read_frame(frame).particles.position[:, 2],
        'moment_x': traj.read_frame(frame).particles.moment_inertia[:, 0],
        'moment_y': traj.read_frame(frame).particles.moment_inertia[:, 1],
        'moment_z': traj.read_frame(frame).particles.moment_inertia[:, 2],
        'angmom_r': traj.read_frame(frame).particles.angmom[:, 0],
        'angmom_i': traj.read_frame(frame).particles.angmom[:, 1],
        'angmom_j': traj.read_frame(frame).particles.angmom[:, 2],
        'angmom_k': traj.read_frame(frame).particles.angmom[:, 3],
        'ori_r': traj.read_frame(frame).particles.orientation[:, 0],
        'ori_i': traj.read_frame(frame).particles.orientation[:, 1],
        'ori_j': traj.read_frame(frame).particles.orientation[:, 2],
        'ori_k': traj.read_frame(frame).particles.orientation[:, 3],
        'v_x': traj.read_frame(frame).particles.velocity[:, 0],
        'v_y': traj.read_frame(frame).particles.velocity[:, 1],
        'v_z': traj.read_frame(frame).particles.velocity[:, 2],
    }
    return pd.DataFrame(particle_info)


def get_bonds(traj, frame):
    bond_info = {
        'type': [traj.read_frame(frame).bonds.types[typeid] 
                 for typeid in traj.read_frame(frame).bonds.typeid],
        'bond_0': traj.read_frame(frame).bonds.group[:, 0],
        'bond_1': traj.read_frame(frame).bonds.group[:, 1]
    }
    return pd.DataFrame(bond_info)


try:
    config = ConfigDict(sys.argv[1])  # Read input parameter file
except:
    print("Usage: %s <config_file> <input_gsd>" % sys.argv[0])
    raise

f = gsd.hoomd.open(sys.argv[2])
PS_spacing = int(config['Packaging Site Spacing'])
all_particles = get_particles(f, 0)
all_bonds = get_bonds(f, 0)

particle_types = np.array(all_particles['type'])
P_types = particle_types[np.where(particle_types == 'P')]
non_P_types = particle_types[np.where(particle_types != 'P')]
P_types[0:len(P_types):PS_spacing] = 'PS'
updated_particle_types = np.concatenate((P_types, non_P_types))
all_particles['type'] = updated_particle_types

# Add packaging site receptors at the base of each ARM
# (Overlapping the single R per ARM that belongs to the rigid body)
PSR_clones = pd.DataFrame(all_particles[all_particles['body']!=-1][all_particles['type']=='R'])
PSR_clones['type'] = ['PSR'] * len(PSR_clones['type'])
all_particles = pd.concat((all_particles, PSR_clones), ignore_index=True)

particle_types = np.unique(all_particles['type']).tolist()
all_particles['typeid'] = [particle_types.index(t) for t in all_particles['type']]
bond_types = np.unique(all_bonds['type']).tolist()
all_bonds['typeid'] = [bond_types.index(t) for t in all_bonds['type']]

box = np.array([40, 40, 40])

# Create HOOMD GSD Snapshot
s = gsd.hoomd.Snapshot()
s.configuration.box = [box[0], box[1], box[2], 0, 0, 0]

# Particle information
s.particles.N = len(all_particles)
s.particles.types = particle_types
s.particles.position = np.array(
    [all_particles['position_x'].tolist(), 
     all_particles['position_y'].tolist(), 
     all_particles['position_z'].tolist()]).T
s.particles.orientation = np.array(
    [all_particles['ori_r'].tolist(), 
     all_particles['ori_i'].tolist(), 
     all_particles['ori_j'].tolist(),
     all_particles['ori_k'].tolist()]).T
s.particles.angmom = np.array(
    [all_particles['angmom_r'].tolist(), 
     all_particles['angmom_i'].tolist(), 
     all_particles['angmom_j'].tolist(),
     all_particles['angmom_k'].tolist()]).T
s.particles.typeid = np.array(all_particles['typeid'])
s.particles.mass = np.array(all_particles['mass'])
s.particles.diameter = np.array(all_particles['diameter'])
s.particles.body = np.array(all_particles['body'])
s.particles.moment_inertia = np.array(
    [all_particles['moment_x'].tolist(), 
     all_particles['moment_y'].tolist(), 
     all_particles['moment_z'].tolist()]).T
s.particles.velocity = np.array(
    [all_particles['v_x'].tolist(), 
     all_particles['v_y'].tolist(), 
     all_particles['v_z'].tolist()]).T

# Bond information
s.bonds.N = len(all_bonds)
s.bonds.types = bond_types
s.bonds.typeid = np.array(all_bonds['typeid'])
s.bonds.group = np.array(
    [all_bonds['bond_0'].tolist(), 
     all_bonds['bond_1'].tolist()]).T

# Write snapshot to initial gsd file
f = gsd.hoomd.open(name='assembled_capsid_PS.gsd', mode='wb')
f.append(s)


