import os
import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from graphenemodeling import graphene

mlg=graphene.Monolayer()
eF = 0.4*constants.e # Coulumbs
kF = mlg.kFermi(eF,model='LowEnergy')
k = np.linspace(-2*kF,2*kF,num=100)

conduction_band = mlg.DiracFermionDispersion(k,model='LowEnergy')
valence_band = -conduction_band

if False:
    fig, ax = plt.subplots(figsize=(5,6))
    ax.plot(k/kF,conduction_band/eF,'k')
    ax.plot(k/kF,valence_band/eF, 'k')
    ax.plot(k/kF,np.zeros_like(k),color='gray')
    ax.axvline(x=0,ymin=0,ymax=1,color='gray')
    ax.set_axis_off()
    plt.show()

if True:
    mlg = graphene.Monolayer()
    mlg.g0prime = -0.2*mlg.g0
    kmax = np.abs(mlg.K)
    emax = mlg.DiracFermionDispersion(0,model='FullTightBinding')
    kx = np.linspace(-kmax,kmax,num=100)
    ky = np.copy(kx)

    # k is relative to K. Add K to move to center of Brillouin zone
    k = (kx + 1j*ky[:,np.newaxis]) + mlg.K

    conduction_band = mlg.DiracFermionDispersion(k,model='FullTightBinding',eh=1)
    valence_band = mlg.DiracFermionDispersion(k,model='FullTightBinding',eh=-1)

    fig = plt.figure(figsize=(8,8))
    fullax = plt.axes(projection='3d')
    fullax.view_init(20,25)
    KX, KY = np.meshgrid(kx,ky)
    fullax.plot_surface(KX/kmax,KY/kmax,conduction_band/mlg.g0,
                rstride=1,cstride=1,cmap='viridis',edgecolor='none')
    fullax.plot_surface(KX/kmax,KY/kmax,valence_band/mlg.g0,
                rstride=1,cstride=1,cmap='viridis',edgecolor='none')

    fullax.set_xlabel('$k_x/|K|$')
    fullax.set_ylabel('$k_y/|K|$')
    fullax.set_zlabel('$\epsilon/\gamma_0$')
    fullax.set_title('Brillouin Zone of Graphene')

    dirname=os.path.dirname(os.path.dirname(__file__))
    plt.savefig(os.path.join('notebooks','Graphene','images','overview.png'))
    plt.show()
