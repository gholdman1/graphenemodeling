"""
================================
Limits (:mod:`photonics.limits`)
================================

"""

def MaterialResponse(OpticalConductivity,omega,d=None,method='svd',restype='material'):
    '''
    A general material optical response limit independent of geometrical parameters.

    See equation 6 of Ref 1.

    Assumes a linear conductivity model K=sigma * E, could be generalized
    to K = L E where L is a differential operator.

    For LDOS options, assumes a planar geometry. See Ref 1 for solid angle correction for nonplanar

    Parameters
    ----------

    OpticalConductivity:

    omega:      array-like, the frequency range of interest (rad/s)

    d:          scalar, Distance from graphene (m), if using LDOS restype (see below)

    method:     'svd': Use singular value decomposition version of formula
                'scalar': Uses the simplified scalar version.

    restype:    Bound on response types given below
                'material': Instrinsic material FOM
                'CSabs': Absorption cross section / area
                'CSext': Extinction cross section / area
                'CSsca': Scattering cross section / area
                'LDOStot': Total LDOS / Free space LDOS
                'LDOSnrad': Nonradiative LDOS / Free space LDOS
                'LDOSrad': Radiative LDOS / Free space LDOS


    References
    ----------

    [1] Miller et al. 2017
        URL: http://pubs.acs.org/doi/10.1021/acs.nanolett.7b02007
    '''

    if np.all(d):
        k0 = omega/sc.speed_of_light
        LDOSprop = 3 / (8 * (k0*d)**4)

    prop = dict({'material':1,
                'CSabs':1, 'CSext': 1, 'CSsca': 1/4,
                'LDOStot':1*LDOSprop,'LDOSnrad':1*LDOSprop,'LDOSrad':LDOSprop/4})

    if method == 'scalar':
        sigma = OpticalConductivity(0,omega,gamma,FermiLevel,T)

        bound = Z_0 * np.abs(sigma)**2 / np.real(sigma)

    elif method == 'svd':
        bound = np.empty_like(omega)

        for i, w in np.ndenumerate(omega):
            s = OpticalConductivityMatrix(0,w,gamma,FermiLevel,T)
            s_dag = np.conj(np.transpose(s))
            s_re_inv = np.linalg.inv(np.real(s))
            prod = np.dot(s_dag,np.dot(s_re_inv,s))
            two_norm= np.linalg.svd(prod)[1][0]
            bound[i] = Z_0*two_norm



    return prop[restype] * bound

