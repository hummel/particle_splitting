# split_part.py
# Jacob Hummel
# created 4/30/15
# Particle splitting routine for gadget files in the hdf5 format.  

import os
import numpy as np
import h5py



def select_region(low_res, header, rmax, nnew):
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
    #r = np.sqrt(x*x + y*y + z*z) / hubble
    #to_refine = np.where(r <= rsplit)[0] #indices of particles to refine
    to_refine = np.where((np.abs(x) <= rmax/2) &
                         (np.abs(y) <= rmax/2) &
                         (np.abs(z) <= rmax/2))[0]
    print to_refine.size, 'particles to refine.'
    return to_refine

def distribute_positions(low_res, hi_res, refine_idx, nnew):
    '''
    Randomly distribute split particles within smoothing kernel
    of unsplit particle, then adjust smoothing length accordingly.
    '''
    pos = low_res['Coordinates'].value
    # Duplicate base position for randomizing
    pos_new = np.repeat(pos[refine_idx], nnew-1, axis=0)
    # randomize position within smoothing kernel of each particle
    hsml = low_res['SmoothingLength'].value
    delta = np.random.rand(refine_idx.size*(nnew-1), 3)
    delta *=2; delta -=1 #convert random sample from (0,1) to (-1,1)
    delta = delta * np.repeat(hsml[refine_idx], nnew-1)[:, np.newaxis]
    pos_new = pos_new + delta
    # Now, randomize original particle position.
    delta = np.random.rand(refine_idx.size, 3)
    delta *=2; delta -=1 #convert random sample from (0,1) to (-1,1)
    delta = delta * hsml[refine_idx, np.newaxis]
    pos[refine_idx] += delta
    # Combine and save.
    pos = np.concatenate((pos, pos_new))
    print pos.shape[0], 'particle positions.'
    hi_res.create_dataset('Coordinates', data=pos)

    # Now shrink smoothing lengths accordingly:
    hsml[refine_idx] *= nnew**(-1./3.)
    hsml = np.concatenate((hsml, np.repeat(hsml[refine_idx], nnew-1)))
    print hsml.size, 'smoothing length entries.'
    hi_res.create_dataset('SmoothingLength', data=hsml)


def refine_mass(low_res, hi_res, refine_idx, nnew):
    '''
    Spread mass evenly over split particles.
    '''
    mass = low_res['Masses'].value
    print 'Mass before refinement:', mass.sum()
    mass[refine_idx] /= nnew
    mass = np.concatenate((mass, np.repeat(mass[refine_idx], nnew-1)))
    print mass.size, 'particle masses.'
    print 'Mass after refinement:', mass.sum()
    hi_res.create_dataset('Masses', data=mass)


def extend_particleIDs(low_res, hi_res, id_max, npart_new):
    pid = low_res['ParticleIDs'].value
    old_max = np.max((pid.max(), id_max))
    id_new = np.arange(old_max+1, old_max+npart_new+1, dtype=pid.dtype)
    pid = np.concatenate((pid,id_new))
    print pid.size, 'particle IDs;', id_new.size, 'new entries.'
    hi_res.create_dataset('ParticleIDs', data=pid)


def extend_field(field, low_res, hi_res, refine_idx, nnew):
    '''
    For all other particle properties, just replicate.
    '''
    values = low_res[field].value
    values = np.concatenate((values, np.repeat(values[refine_idx],
                                               nnew-1, axis=0)))
    hi_res.create_dataset(field, data=values)


def split_particles(low_res, hi_res, header, id_max, rsplit, nnew,):
    refine_idx = select_region(low_res, header, rsplit, nnew)
    # Spread particles over smoothing kernel.
    distribute_positions(low_res, hi_res, refine_idx, nnew)
        # Spread mass evenly over split particles.
    refine_mass(low_res, hi_res, refine_idx, nnew)
    # Extend ParticleIDs to cover new particles.
    nnew_ids = refine_idx.size*(nnew-1)
    extend_particleIDs(low_res, hi_res, id_max, nnew_ids)
    # Extend all other particle fields.
    all_fields = low_res.keys()
    special = ['ParticleIDs', 'Masses', 'Coordinates', 'SmoothingLength']
    others = [x for x in all_fields if x not in special]
    for field in others:
        extend_field(field, low_res, hi_res, refine_idx, nnew)
    #Update particle counts in header.
    npart_total = header.attrs['NumPart_Total']
    npart_this_file = header.attrs['NumPart_ThisFile']
    npart_total[0] += nnew_ids
    npart_this_file[0] += nnew_ids
    header.attrs.modify('NumPart_Total', npart_total)
    header.attrs.modify('NumPart_ThisFile', npart_this_file)


def main(filein, fileout, rmax, nnew, includeDM=True, physical_units=False):
    print "This is split_part.py."
    print "Splitting gas particles into", nnew, "daughters."
    infile = h5py.File(filein, 'r')
    outfile = h5py.File(fileout, 'w') 
    # Duplicate header from infile to outfile.
    header = outfile.create_group('Header')
    for entry in infile['Header'].attrs.items():
        header.attrs.create(*entry)

    # Convert from physical units to comoving if desired.
    if physical_units:
        print "Splitting inner", rmax, "physical kpc."
        h = header.attrs['HubbleParam']
        a = header.attrs['Time']
        rmax = rmax * h / a
    print "Splitting inner", rmax, "comoving kpc/h of", filein

    # Don't want to split dark matter particles, so simply duplicate the data.
    infile.copy('PartType1', outfile)
    dmid_max = infile['PartType1/ParticleIDs'].value.max()

    # Split gas particles
    old_gas = infile['PartType0']
    new_gas = outfile.create_group('PartType0')
    new_gas = split_particles(old_gas, new_gas, header, dmid_max, rmax, nnew)

    infile.close()
    outfile.close()
    print '\n\n'

if __name__ == '__main__':
    infile = os.getenv("HOME")+"/sim/stampede/halo2_vanilla/snapshot_0145.hdf5"
    outfile = os.getenv("HOME")+"/sim/stampede/halo2_vanilla/snapshot_0146.hdf5"
    outfile2 = os.getenv("HOME")+"/sim/stampede/halo2_vanilla/snapshot_0147.hdf5"
    #rmax = 10
    #nnew = 8
    main(infile, outfile, 8,8, physical_units=False)
    main(outfile, outfile2,6,8, physical_units=False)
