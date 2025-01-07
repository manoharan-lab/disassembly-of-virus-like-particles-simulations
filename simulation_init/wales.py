import os
import sys
import math
from config_dict import ConfigDict
import hoomd
import hoomd.md as md

try:
    config = ConfigDict(sys.argv[1])  # Read input parameter file
    seed = int(sys.argv[2])           # Get unique random seed
except:
    print("Usage: %s <input_file> <seed>" % sys.argv[0])
    raise

# Load initial state
hoomd.context.initialize("")

if os.path.isfile('traj.gsd'):
    hoomd.init.read_gsd('traj.gsd', frame=-1)
else:
    hoomd.init.read_gsd('assembled_capsid_PS.gsd', time_step=0)

# Define capsomer rigid bodies
rigid = hoomd.md.constrain.rigid()
rigid.set_param('M',
                types=['C', 'C', 'C', 'C', 'C',
                       'A', 'A', 'A', 'A', 'A',
                       'T', 'B',
                       'X', 'X', 'X', 'X', 'X',
                       'X', 'X', 'X', 'X', 'X',
                       'X', 'R', 'X', 'R', 'X',
                       'R', 'X', 'R', 'X', 'R',
                       'X', 'X', 'X', 'X', 'X',
                       'PSR', 'PSR', 'PSR', 'PSR', 'PSR'],
                positions=[(0.375447, -0.530596, -0.255473),
                           (-0.388553, -0.521073, -0.255393),
                           (-0.615601, 0.208651, -0.255357),
                           (0.00804895, 0.649958, -0.255396),
                           (0.620652, 0.193124, -0.25535),
                           (0.57765, -0.81633, 0.0945986),
                           (-0.597865, -0.801604, 0.0945797),
                           (-0.947164, 0.320865, 0.0946718),
                           (0.0124426, 0.999954, 0.0946254),
                           (0.95482, 0.297107, 0.0946505),
                           (0.0000146141, -0.0000253554, 0.594546),
                           (0.0000146141, -0.0000253554, -0.405433),
                           (0.481359, -0.680251, 0.094512),
                           (-0.0084051, -0.674148, 0.0945124),
                           (-0.498179, -0.668035, 0.094518),
                           (-0.643747, -0.200365, 0.0945791),
                           (-0.789227, 0.267415, 0.0946146),
                           (-0.389457, 0.55027, 0.0947003),
                           (0.0103689, 0.833226, 0.0946563),
                           (0.403015, 0.540477, 0.094611),
                           (0.795701, 0.247638, 0.0946585),
                           (0.638532, -0.216351, 0.0945892),
                           (0.288869, -0.408183, 0.0945285),
                           (0.288776, -0.408108, -0.255456),
                           (-0.298937, -0.400809, 0.0945462),
                           (-0.298903, -0.400846, -0.255394),
                           (-0.473547, 0.160415, 0.0945709),
                           (-0.4735, 0.160489, -0.255365),
                           (0.00619251, 0.499969, 0.094628),
                           (0.00619111, 0.499943, -0.255396),
                           (0.477409, 0.148591, 0.0945787),
                           (0.477444, 0.148565, -0.25536),
                           (0.0962862, -0.136005, 0.0945854),
                           (-0.0996033, -0.133583, 0.094536),
                           (-0.157804, 0.0535246, 0.0946172),
                           (0.00210007, 0.166603, 0.0945666),
                           (0.159159, 0.0495531, 0.0945922),
                           (0.288776, -0.408108, -0.255456),
                           (-0.298903, -0.400846, -0.255394),
                           (-0.4735, 0.160489, -0.255365),
                           (0.00619111, 0.499943, -0.255396),
                           (0.477444, 0.148565, -0.25536)]
                )
rigid.validate_bodies()

# Get simulation parameters
rho = float(config['rho'])
subunit_attraction = float(config['Ebond'])

# Polymer (P) interactions
pol_sigma = 0.1
pol_bond_length = pol_sigma
pol_cut = pol_sigma
pol_bend_strength = float(config['Pol Bend Strength'])  # spring coefficient of bending
packaging_site_attraction = float(config['Packaging Site Strength'])

# Top (T) interactions
top_sigma = 2.1
top_cut = top_sigma
#top_repulsion = subunit_attraction / 4.0
top_repulsion = 3.0 / 4.0  # Subunit attraction for assembly

# Bottom (B) interactions
bot_sigma = 1.8
bot_cut = bot_sigma

# Excluder (X) interactions
x_sigma = 0.35
x_cut = x_sigma

unit_length = 5  # 1 length unit = 5nm
debye_length = float(config['Debye'])
kappa = 1 / debye_length
bjerrum_length = 0.7 / unit_length  # 0.7nm in water at room temperature
yukawa_prefactor = bjerrum_length \
                   * debye_length \
                   * math.exp(pol_sigma / debye_length) / (pol_sigma + debye_length)  # magnitude of yukawa potential
# set a minimum cutoff radius
yukawa_rcut = 3 * max(debye_length, 0.2)   # Cutoff radius beyond which Yukawa potential is zero
yukawa_ron = 2 * max(debye_length, 0.2)    # Onset of smoothing to zero

'''
Pseudoatom types

A: Attractor (subunit-subunit attraction)
T: Top of capsomer (repulsive)
B: Bottom of capsomer (repulsive with top)
P: Polymer
PS: Packaging signal (there are none)
R: ARM charged component
RA: (not used)
N: ARM neutral component
Q: (there are none)
X: Excluder
C: Something that only interacts with PS's

'''

nl = md.nlist.cell(check_period=1)  # Cell-based neighbor list

# Lennard-Jones potential: 4 * epsilon [(sigma/r)**12 - alpha * (sigma/r)**6]
lj = md.pair.lj(r_cut=top_cut, nlist=nl)
lj.set_params(mode="shift")  # A constant shift is applied to the entire potential so that it is 0 at the cutoff

lj.pair_coeff.set('T', 'T',
                  alpha=0, epsilon=top_repulsion, r_cut=top_cut, sigma=top_sigma)
lj.pair_coeff.set('B', 'T',
                  alpha=0, epsilon=top_repulsion, r_cut=bot_cut, sigma=bot_sigma)
lj.pair_coeff.set(['P', 'PS', 'R', 'RA', 'N', 'Q'], ['P', 'PS', 'R', 'RA', 'Q'],
                  alpha=0, epsilon=1, r_cut=pol_cut, sigma=pol_sigma)
lj.pair_coeff.set(['P', 'PS', 'R', 'RA', 'N', 'Q', 'X'], 'X',
                  alpha=0, epsilon=1, r_cut=x_cut, sigma=x_sigma)
lj.pair_coeff.set(['M', 'P', 'PS', 'R', 'RA', 'N', 'Q', 'X', 'C', 'A', 'PSR'],
                  ['M', 'T', 'B', 'N', 'A', 'C', 'PSR'],
                  alpha=0, epsilon=0, r_cut=0, sigma=0)
lj.pair_coeff.set('B', 'B',
                  alpha=0, epsilon=0, r_cut=0, sigma=0)

# Yukawa potential: epsilon * exp(-kappa * r) / r
yukawa = md.pair.yukawa(r_cut=yukawa_rcut, nlist=nl)
yukawa.set_params(mode="xplor")  # use a smoothing function between r_on and r_cut

yukawa.pair_coeff.set(['P', 'PS'], ['R', 'RA'],
                      epsilon=-yukawa_prefactor, kappa=kappa, r_cut=yukawa_rcut, r_on=yukawa_ron)
yukawa.pair_coeff.set(['P', 'PS'], ['P', 'PS'],
                      epsilon=yukawa_prefactor, kappa=kappa, r_cut=yukawa_rcut, r_on=yukawa_ron)
yukawa.pair_coeff.set(['R', 'RA'], ['R', 'RA'],
                      epsilon=yukawa_prefactor, kappa=kappa, r_cut=yukawa_rcut, r_on=yukawa_ron)
yukawa.pair_coeff.set(['M', 'C', 'A', 'B', 'T', 'X', 'P', 'PS', 'R', 'RA', 'N', 'Q', 'X3', 'PSR'],
                      ['M', 'C', 'A', 'B', 'T', 'X', 'X3', 'N', 'PSR'],
                      epsilon=0, kappa=0, r_cut=yukawa_rcut, r_on=yukawa_ron)

# Morse potential: D0[exp(-2 alpha(r-r0)) - 2 exp(-alpha(r-r0))]
r0 = 0.2
alpha = rho / r0
morse = md.pair.morse(r_cut=2.0, nlist=nl)
morse.set_params(mode="shift")
morse.pair_coeff.set(['M', 'C', 'A', 'B', 'T', 'X', 'P', 'PS', 'R', 'RA', 'N', 'Q', 'X3', 'PSR'],
                     ['M', 'C', 'B', 'T', 'X', 'X3', 'P', 'R', 'RA', 'N', 'Q'],
                     D0=0, alpha=0, r0=0, r_cut=0.0)
morse.pair_coeff.set('A', 'A',
                     D0=subunit_attraction, alpha=alpha, r0=r0, r_cut=2.0)
morse.pair_coeff.set(['PSR', 'PS'], 'A',
                     D0=0, alpha=0, r0=0, r_cut=0.0)
morse.pair_coeff.set('PSR', 'PSR',
                     D0=0, alpha=0, r0=0, r_cut=0.0)
morse.pair_coeff.set('PS', 'PS',
                     D0=0, alpha=0, r0=0, r_cut=0.0)
morse.pair_coeff.set('PSR', 'PS',
                     D0=packaging_site_attraction, alpha=alpha, r0=r0, r_cut=2.0)

# Polymer bonds
harmonic = md.bond.harmonic()
harmonic.bond_coeff.set('polymer', k=330.0 * unit_length**2, r0=pol_bond_length)
harmonic.bond_coeff.set('tail1', k=330.0 * unit_length**2, r0=0.35)
harmonic.bond_coeff.set('tail2', k=330.0 * unit_length**2, r0=pol_bond_length)
if pol_bend_strength > 0.0:
    harmonic = md.angle.harmonic()
    harmonic.angle_coeff.set('polymer', k=pol_bend_strength, t0=math.pi)

nl.reset_exclusions(exclusions=['bond', 'body'])

# Set the integrator
ts = 0.0001
animation_period = int(2.0 / ts)

group_rigid = hoomd.group.rigid_center()
group_nonrigid = hoomd.group.nonrigid()
group_integrate = hoomd.group.union('integrate', group_rigid, group_nonrigid)

hoomd.dump.gsd(filename="traj.gsd", period=animation_period, group=hoomd.group.all(), phase=0)
hoomd.dump.dcd(filename="traj.dcd", period=animation_period, group=hoomd.group.all(), phase=0)
hoomd.analyze.log(filename='out.log', quantities=['temperature', 'potential_energy'],
                  period=animation_period, overwrite=False, phase=0)

integrator_mode = md.integrate.mode_standard(dt=ts)
md.integrate.langevin(group=group_integrate, kT=1.0, seed=seed, dscale=True)

num_steps = int(5000.0 / ts)
hoomd.run_upto(num_steps)
