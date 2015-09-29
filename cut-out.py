# cut-out.py
# Jacob Hummel
# created 5/7/15
# Particle splitting / cut-out routine for gadget files in the hdf5 format.  

import os
import numpy as np
import h5py



def select_region(low_res, header, rmax, center=None):
    pos = low_res['Coordinates']
    time = header.attrs['Time']
    hubble = header.attrs['HubbleParam']
    # Center refinement on highest density particle
    if center is None:
        dens = low_res['Density']
        center = pos[dens[:].argmax()]
    x = pos[:,0] - center[0]
    y = pos[:,1] - center[1]
    z = pos[:,2] - center[2]
    # Calculate distance from center in physical parsecs
    #r = np.sqrt(x*x + y*y + z*z) / hubble
    #to_refine = np.where(r <= rmax)[0] #indices of particles to refine
    to_refine = np.where((np.abs(x) <= rmax/2) &
                         (np.abs(y) <= rmax/2) &
                         (np.abs(z) <= rmax/2))[0]
    print to_refine.size, 'particles to refine.'
    return to_refine, center

def distribute_gas(low_res, hi_res, refine_idx, nnew, center, rmax):
    '''
    Randomly distribute split particles within smoothing kernel
    of unsplit particle, then adjust smoothing length accordingly.
    '''
    pos = low_res['Coordinates'].value - center + rmax/2
    # Duplicate base position for randomizing
    pos= np.repeat(pos[refine_idx], nnew, axis=0)
    # randomize position within smoothing kernel of each particle
    hsml = low_res['SmoothingLength'].value
    delta = np.random.rand(refine_idx.size*nnew, 3)
    delta *=2; delta -=1 #convert random sample from (0,1) to (-1,1)
    delta = delta * np.repeat(hsml[refine_idx], nnew)[:, np.newaxis]
    pos += delta
    print pos.shape[0], 'particle positions.'
    #pos -= pos.min()
    hi_res.create_dataset('Coordinates', data=pos)

    # Now shrink smoothing lengths accordingly:
    hsml[refine_idx] *= nnew**(-1./3.)
    hsml = np.repeat(hsml[refine_idx], nnew)
    print hsml.size, 'smoothing length entries.'
    hi_res.create_dataset('SmoothingLength', data=hsml)


def refine_mass(low_res, hi_res, refine_idx, nnew):
    '''
    Spread mass evenly over split particles.
    '''
    mass = low_res['Masses'].value
    print 'Total mass before refinement:', mass.sum()
    print 'Particle mass before refinement:', mass.min()
    mass[refine_idx] /= nnew
    mass = np.repeat(mass[refine_idx], nnew)
    print mass.size, 'particle masses.'
    print 'Total mass after refinement:', mass.sum()
    print 'Particle mass after refinement:', mass.min()
    hi_res.create_dataset('Masses', data=mass)


def extend_particleIDs(low_res, hi_res, id_max, npart_new):
    pid = low_res['ParticleIDs'].value
    id_new = np.arange(1, npart_new+1, dtype=pid.dtype)
    hi_res.create_dataset('ParticleIDs', data=id_new)


def extend_field(field, low_res, hi_res, refine_idx, nnew):
    '''
    For all other particle properties, just replicate.
    '''
    values = low_res[field].value
    values = np.repeat(values[refine_idx], nnew, axis=0)
    hi_res.create_dataset(field, data=values)


def split_gas(low_res, hi_res, header, id_max, rmax, nnew):
    refine_idx, center = select_region(low_res, header, rmax)
    # Spread mass evenly over split particles.
    refine_mass(low_res, hi_res, refine_idx, nnew)
    # Spread particles over smoothing kernel.
    distribute_gas(low_res, hi_res, refine_idx, nnew, center, rmax)
    # Extend ParticleIDs to cover new particles.
    nnew_ids = refine_idx.size*nnew
    extend_particleIDs(low_res, hi_res, id_max, nnew_ids)
    # Extend all other particle fields.
    all_fields = low_res.keys()
    special = ['ParticleIDs', 'Masses', 'Coordinates', 'SmoothingLength']
    others = [x for x in all_fields if x not in special]
    for field in others:
        extend_field(field, low_res, hi_res, refine_idx, nnew)
    return refine_idx.size * nnew, center

def slice_dm(low_res, hi_res, header, id_max, rmax, center):
    refine_idx, center = select_region(low_res, header, rmax, center)
    n_dm = refine_idx.size
    fields = [x for x in low_res.keys() if x not in ('Coordinates','ParticleIDs')]
    for field in fields:
        values = low_res[field].value
        hi_res.create_dataset(field, data=values[refine_idx])
    pos = low_res['Coordinates'].value - center + rmax/2
    hi_res.create_dataset('Coordinates', data=pos[refine_idx])
    pid = low_res['ParticleIDs'].value
    id_new = np.arange(id_max+1, id_max+n_dm+1, dtype=pid.dtype)
    hi_res.create_dataset('ParticleIDs', data=id_new)
    return n_dm

def main(filein, fileout, rmax, nnew, includeDM=True, physical_units=False):
    print "This is cut-out.py."
    print "Splitting gas particles into", nnew, "daughters."
    infile = h5py.File(filein, 'r')
    outfile = h5py.File(fileout, 'w') 
    # Duplicate header from infile to outfile.
    header = outfile.create_group('Header')
    for entry in infile['Header'].attrs.items():
        header.attrs.create(*entry)

    # Convert from physical units to comoving if desired.
    if physical_units:
        print "Cutting out inner", rmax, "physical kpc."
        h = header.attrs['HubbleParam']
        a = header.attrs['Time']
        rmax = rmax * h / a
    print "Cutting out inner", rmax, "comoving kpc/h of", filein
    # Cut out and split gas particles
    old_gas = infile['PartType0']
    new_gas = outfile.create_group('PartType0')
    npart_gas, center = split_gas(old_gas, new_gas, header, 0, rmax, nnew)
    print "New NumPart_gas:", npart_gas

    # Dark matter cut-out
    new_dm = outfile.create_group('PartType1')
    if includeDM:
        old_dm = infile['PartType1']
        npart_dm = slice_dm(old_dm, new_dm, header, npart_gas, rmax, center)
    else:
        npart_dm = 0
    print "New NumPart_dm:", npart_dm

    #Update particle counts in header.
    print "New NumPart_Total:", npart_gas + npart_dm
    npart_total = header.attrs['NumPart_Total']
    npart_this_file = header.attrs['NumPart_ThisFile']
    npart_total[0] = npart_gas
    npart_this_file[0] += npart_gas
    npart_total[1] = npart_dm
    npart_this_file[1] += npart_dm
    header.attrs.modify('NumPart_Total', npart_total)
    header.attrs.modify('NumPart_ThisFile', npart_total)


    infile.close()
    outfile.close()
    print '\n\n'

if __name__ == '__main__':
    infile = os.getenv("HOME")+"/sim/lonestar/halo2_lowres/snapshot_0173.hdf5"
    outfile = os.getenv("HOME")+"/sim/halo2/small/snapshot_0174.hdf5"
    main(infile, outfile, .01,8, includeDM=True, physical_units=True)
