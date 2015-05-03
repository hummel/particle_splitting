# split_part.py
# Jacob Hummel
# created 4/30/15
# Particle splitting routine for gadget files in the hdf5 format.  

import os
import numpy as np
import h5py



def select_region(low_res, header, rsplit, nnew):
	pos = low_res['Coordinates']
	dens = low_res['Density']
	time = header.attrs['Time']
	hubble = header.attrs['HubbleParam']
	# Center refinement on highest density particle
	center = pos[dens[:].argmax()]
	x = pos[:,0] - center[0]
	y = pos[:,1] - center[1]
	z = pos[:,2] - center[2]
	# Calculate distance from center in physical parsecs
	r = np.sqrt(x*x + y*y + z*z) / hubble
	to_refine = np.where(r <= rsplit)[0] #indices of particles to refine
	print to_refine.size, 'particles to refine.'
	return to_refine

def refine_mass(low_res, hi_res, refine_idx):
	# Spread mass evenly over split particles.
	mass = low_res['Masses'].value
	for ix in refine_idx:
		mass[ix] /= nnew
	mass = np.concatenate((mass, np.repeat(mass[refine_idx], nnew-1)))
	hi_res.create_dataset('Masses', data=mass)

def extend_field(field, low_res, hi_res, refine_idx):
	# For all other particle properties, just replicate.
	 values = low_res[field].value
	 if values.ndim == 1:
		 values = np.concatenate((values, np.repeat(values[refine_idx], nnew-1)))
		 hi_res.create_dataset(field, data=values)

def split_particles(low_res, hi_res, header, rsplit, nnew):
	refine_idx = select_region(low_res, header, rsplit, nnew)
	pos = low_res['Coordinates'].value
	print pos.max(), pos.min()
	print pos[refine_idx].max(), pos[refine_idx].min()	
	refine_mass(low_res, hi_res, refine_idx)

	all_fields = low_res.keys()
	special = ['Masses', 'Coordinates', 'ParticleIDs', 'SmoothingLength']
	others = [x for x in all_fields if x not in special]
	for field in others:
		extend_field(field, low_res, hi_res, refine_idx)


def main(infile, outfile, rsplit, nnew):
	infile = h5py.File(infile, 'r')
	outfile = h5py.File(outfile, 'w') 
	# Duplicate header from infile to outfile.
	header = outfile.create_group('Header')
	for entry in infile['Header'].attrs.items():
		header.attrs.create(*entry)
	# Don't want to split dark matter particles, so simply duplicate the data.
	infile.copy('PartType1', outfile)
	
	# Split gas particles
	old_gas = infile['PartType0']
	new_gas = outfile.create_group('PartType0')
	new_gas = split_particles(old_gas, new_gas, header, rsplit, nnew)
	
	infile.close()
	outfile.close()

if __name__ == '__main__':
	infile = os.getenv("HOME")+"/sim/lonestar/vanilla2/snapshot_0176.hdf5"
	outfile = os.getenv("HOME")+"/sim/duplicate.hdf5"
	rsplit = 10
	nnew = 8
	main(infile, outfile, rsplit, nnew)
