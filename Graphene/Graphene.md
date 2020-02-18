*See zotero for latest version*

# Graphene

[TOC]

# Band Structure

## Linear Dispersion

Valid for wavevectors satisfying $\hbar v_F k_F \ll 8.25$ eV, since that is the high-energy cutoff of the lattice (i.e. the energy at the edge of the Brillouin Zone).

# Density-Density Response

Many of the properties of graphene can be recovered from the density-density response $\chi_0(\mathbf{r,r'};\omega)$.

## Conductivity

*See Graphene-OpticalProperties*

The optical conductivity of monolayer graphene can be modeled with the [RPA](https://en.wikipedia.org/wiki/Random_phase_approximation) as:
$$
\sigma_s(\omega) = \frac{i2e^2k_BT}{\pi\hbar^2(\omega+i\tau^{-1})}\ln\left[ 2 \cosh\left(\frac{E_F}{2k_BT}\right)  \right]\\ + \frac{e^2}{4\hbar}\left[ \frac{1}{2} + \frac{1}{\pi}\arctan\left(\frac{\hbar\omega - 2E_F}{2k_BT}\right)\\ - \frac{i}{2\pi}\ln\frac{(\hbar\omega+2E_F)^2}{(\hbar\omega-2E_F)^2+4(k_BT)^2} \right]
$$
Originally derived in [Falkovsky and Pershoguba 2007](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.76.153410) and in used in

- [Yao et al. 2013](https://pubs.acs.org/doi/10.1021/nl3047943)
- [Li and Yu 2013](https://aip.scitation.org/doi/10.1063/1.4800931)

The [Photonics 101](http://photonics101.com/light-matter-interactions/the-permittivity-of-graphene) website gives
$$
\sigma_{2D}(\omega) = i \frac{1}{\pi\hbar^2}\frac{e^2k_BT}{\omega + i2\Gamma}\left\{ \frac{\mu_c}{k_BT} + 2\ln \left[ \exp(-\mu_c/k_BT) + 1\right] \right\}\\ + i \frac{e^2}{4\pi\hbar}\ln\left[ \frac{2|\mu_c|-\hbar(\omega + i2\Gamma)}{2|\mu_c|+\hbar(\omega+i2\Gamma)} \right]
$$
which may be the same, but I have yet to compare. 

Note: There seems to a factor of 2 difference between $\tau^{-1}$ and $\Gamma$. I believe this may be related to the chosen definition for each. In one case, the lifetime may refer to the lifetime of the amplitude A of an oscillation, while in another case it may refer to the lifetime of the intensity $\propto |A|^2$.

Other places to find conductivity of graphene

* [Self-consistent field approach](https://link.aps.org/doi/10.1103/PhysRev.115.786)
* [Kubo formalism](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.266802)
* RPA: PM Platzman and PA Wolff *Waves and Interactions in Solid-State Plasma*

### Drude/Intraband Component

 From [Gusynin, Sharapov, and Carbotte 2006](http://iopscience.iop.org/article/10.1088/0953-8984/19/2/026222/meta)
$$
\sigma_{xx}^{\text{Drude}}(\Omega) = -\frac{2ie^2(\Omega+2i\Gamma)}{h}\frac{1}{(\Omega+2i\Gamma)^2}\int_\Delta^\infty d\omega \frac{\omega^2-\Delta^2}{\omega}\left( \frac{\partial n_F(\omega)}{\partial\omega} - \frac{\partial n_F(\omega)}{\partial\omega} \right)
$$

where $\Delta$ gives the asymmetry between sublattices. Lets assume for now that $\Delta=0$. To evaluate the integral, we note that ...



### Full Derivation

#### From Kubo Formula

Originally derived by [Gusynin, Sharapov, and Carbotte 2006](http://iopscience.iop.org/article/10.1088/0953-8984/19/2/026222/meta). Notes given here are meant to supplement the derivation in the paper. Much of the phrasing is directly quoted from the paper since I am unfamiliar with many of the terms. Watch out for missing factors of $\hbar$ or $k_B$.

The Kubo Formula is given by
$$
\sigma_{ij} = \frac{\prod _{ij}^R(\Omega + i0)}{i\Omega}
$$
where the retarded current-current correlation function in the bubble approximation is given by
$$
\prod_{ij}(\Omega+i0) = e^2v_F^2\int_{-\infty}^\infty d\omega d\omega' \frac{n_F(\omega')-n_F(\omega)}{\omega-\omega'-\Omega-i0} \times \int \frac{d^2k}{(2\pi)^2}\text{tr}\left[ \gamma^iA(\omega,\mathbf{k})\gamma^jA(\omega',\mathbf{k}) \right]
$$
where $n_F(\omega)= [\exp((\hbar\omega-\mu)/k_BT)+1]^{-1}$ is the Fermi function and
$$
\text{Spectral Function: }A(\omega,\mathbf{k})=e^{-ck^2/|eB|}\sum_{n=0}^\infty \frac{(-1)^n\Gamma_n(\omega)}{2\pi M_n}\left[ \frac{(\gamma^0M_n+\Delta)f_1(\mathbf{k})+f_2(\mathbf{k})}{(\omega-M_n)^2+\Gamma_n^2(\omega)} +\\
\frac{(\gamma^0M_n-\Delta)f_1(\mathbf{k})-f_2(\mathbf{k})}{(\omega+M_n)^2+\Gamma_n^2(\omega)}\right]
$$

$$
\text{Index of Landau Level: } n
$$

$$
f_1(\mathbf{k}) = 2\left[ P_-L_n\left(\frac{2ck^2}{|eB|}\right) - P_+L_{n-1}\left(\frac{2ck^2}{|eB|}\right) \right] \\
f_2(\mathbf{k}) = 4v_F\mathbf{k}\gamma L_{n-1}^1\left(\frac{2ck^2}{|eB|}\right)
$$

$$
\text{Scattering Rate: } \Gamma_n(\omega) \approx \Gamma(\omega)
$$

$$
\text{Projectors: }P_\pm=\frac{1}{2}(1\pm i\gamma^1\gamma^2\text{sign}(eB))
$$

$$
\text{Laguerre Polynomials: } L_n^\alpha(z), L_n(z)=L_n^0(z)
$$

$$
\text{Relativistic Landau Level Energy:} E_n=\pm M_n =\pm\sqrt{\Delta^2+2nv_F^2|eB|/c}
$$

$$
\text{Excitonic Gap: } \Delta = see refs
$$

Here, the spectral function  is "associated with the translation invariant part of the Dirac fermion Green's function in an external magnetic field $\mathbf{B}$ applied perpendicular to the plane along the positive $z$ direction" and is "decomposed over Landau Levels".

##### Low-Field Lorentzian

*Low-Field* is determined by $\sqrt{\hbar|eB|v_F^2/c}<<\Gamma$. Be careful about factors
$$
\boxed{
    \sigma_{xx}(\Omega) = -\frac{2ie^2(\Omega+2i\Gamma)}{h}\left[ \frac{1}{(\Omega+2i\Gamma)^2}\int_\Delta^\infty d\omega \frac{\omega^2-\Delta^2}{\omega^2}\left(\frac{\partial n_F(\omega)}{\partial\omega}-\frac{\partial n_F(-\omega)}{\partial\omega}\right) \right] \\
    -\int_\Delta^\infty d\omega \frac{\omega^2+\Delta^2}{\omega^2}\frac{n_F(-\omega)-n_F(\omega)}{(\Omega+2i\Gamma)^2-4\omega^2}
    \\
    \sigma_{xy}(\Omega) = \frac{e^2v_F^2eB}{\pi c}\int_\Delta^\infty d\omega \left( \frac{\partial n_F(\omega)}{\partial\omega} + \frac{\partial n_F(-\omega)}{\partial\omega}\right)\times \\ \left[ -\frac{\omega^2-\Delta^2}{\omega^2}\frac{1}{(\Omega+2i\Gamma)^2} + \frac{\omega^2+\Delta^2}{\omega^2}\frac{1}{4\omega^2-(\Omega+2i\Gamma)^2} \right]
}
$$

### Approximations

#### Drude

In a Drude model, conductivity is $\sigma(\omega)=\sigma_0/(1+i\omega\tau)​$ where $\tau​$ is the relaxation time and $\sigma_0=n_{3D}e^2\tau/m​$ is the DC conductivity. For frequencies such that $\omega\gg\tau^{-1}​$, this simplifies to $\sigma(\omega)=-ine^2/m\omega​$. For 2D conductivities, we can think of this like $\sigma=-in_{2D}e^2/m\omega​$ The effective mass is $m=p/v_g=\hbar k_F / v_F=\hbar\sqrt{\pi n_{2D}}/v_F​$. So
$$
\sigma(\omega)=-i\frac{e^2v_F\sqrt{n_{2D}}}{\hbar\sqrt{\pi}\omega}
$$

## Landau Levels

Hamiltonian in a magnetic field is given by
$$
H = v_0\sigma\cdot(p-eA)
$$
Leads to levels spaced as
$$
E_n = \pm v_0\sqrt{2eBn}
$$
with eigenstates
$$
\left|n,p_y\right>=e^{-ieBYx}\left(\begin{array}{}\phi^n(y)\\ \pm\phi^{n-1}(y) \end{array}\right)
$$


## Derivation

In Landau gauge, $A=-By\hat x$ so
$$
H = v_0(\sigma_x(p_x+eBy)+\sigma_yp_y)
$$
Define $Y=-\frac{1}{eB}p_x$ to get
$$
H=v_0(\sigma_xBe(y-Y)+\sigma_yp_y)
$$
In this way we can define canonical variables that satisfy commutation relations of harmonic oscillator. The rest follows from there.



## Dipole Polarizability

See equation 5.22 of [this thesis](http://link.springer.com/10.1007/978-3-319-48562-1). The dipole polarizability $\alpha(\omega)​$ is defined by $p(\omega)=\epsilon_0\bar\epsilon\alpha(\omega)E_0​$ and is calculated as 
$$
\alpha(\omega)=2L^3\sum_\nu \frac{|\left<\tilde x|\rho_{||\nu}\right>|^2}{\zeta_\nu-\zeta(\omega)}
$$

See all of chapter 5