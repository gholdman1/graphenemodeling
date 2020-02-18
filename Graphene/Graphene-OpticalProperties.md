# Optical Conductivity

The Optical Conductivity describes the surface current response to an electric field
$$
\mathbf K = \int \sigma(\mathbf r_{||},\mathbf r'_{||};\omega)\mathbf E_{||}(\mathbf r_{||}',\omega)\text{d}^2\mathbf{r_{||}'}
$$
In the local limit,
$$
\mathbf K=\sigma(\omega)\mathbf E_{||}(\omega)
$$


## General Form

Assuming only a tight-binding model, [Falkovsky and Varlamov 2007](https://link.springer.com/article/10.1140/epjb/e2007-00142-3) find
$$
\frac{\sigma_{ij}(\omega,k)}{(ie^2/\pi^2)}=\sum_{a\in bands}\int v^iv^j\frac{f_0(\epsilon_a(\mathbf p_-))-f_0(\epsilon_a(\mathbf p_+))}{(\epsilon_a(\mathbf p_+)-\epsilon_a(\mathbf p_-))(\omega-\epsilon_a(\mathbf p_+)+\epsilon_a(\mathbf p_-))}d^2p\\
+2\omega\int v_{12}^iv_{21}^j\frac{f_0(\epsilon_1(\mathbf p_-))-f_0(\epsilon_2(\mathbf p_+))}{[\epsilon_2(\mathbf p_+)-\epsilon_1(\mathbf p_-)][\omega^2-(\epsilon_a(\mathbf p_+)+\epsilon_a(\mathbf p_-))^2]} d^2p
$$
which does not make any assumptions about graphene except the number of bands. Includes spin and valley symmetry.

## Local Approximation

Discussed in [Christensen 2017](https://link.springer.com/book/10.1007%2F978-3-319-48562-1)

In the limit $k/k_F\to 0​$,
$$
\sigma_{intra}(\omega)=\frac{2ie^2k_BT}{\pi\hbar^2(\omega+i\gamma)}\ln\left[2\cosh\left(\frac{E_F}{2k_BT}\right)\right]
$$

$$
\sigma_{inter}(\omega)=\frac{e^2}{4\hbar}\left[
H(\hbar\omega/2)+\frac{4i}{\pi}\hbar(\omega+i\gamma)\int_{0}^\infty \frac{H(\epsilon)-H(\hbar\omega/2)}{\hbar^2(\omega+i\gamma)^2-4\epsilon^2}d\epsilon
\right]
$$

with 
$$
H(\epsilon)\equiv f(-\epsilon)-f(\epsilon)=\frac{\sinh(\epsilon/k_BT)}{\cosh(E_F/k_BT)+\cosh(\epsilon/k_BT)}
$$

## Local and Low Temperature Approximation

Valid if local and $T<<T_F=E_F/k_B$ (Note: $T_F\approx 4.6\times 10^3 K$ for $E_F=0.4$ eV)
$$
\sigma_{intra}(\omega)=\frac{ie^2E_F}{\pi \hbar^2(\omega+i\gamma)}
$$

$$
\sigma_{inter}(\omega)=\frac{e^2}{4\hbar}\left[\Theta(\hbar\omega-2E_F)+\frac{i}{\pi}\ln\left|\frac{2E_F-\hbar\omega}{2E_F+\hbar\omega}\right| \right]
$$



## Some forms

$$
\sigma_{2D}(\omega) = i \frac{1}{\pi\hbar^2}\frac{e^2k_BT}{\omega + i2\Gamma}\left\{ \frac{\mu_c}{k_BT} + 2\ln \left[ \exp(-\mu_c/k_BT) + 1\right] \right\}\\ + i \frac{e^2}{4\pi\hbar}\ln\left[ \frac{2|\mu_c|-\hbar(\omega + i2\Gamma)}{2|\mu_c|+\hbar(\omega+i2\Gamma)} \right]
$$



# Absorbed Power

The power absorbed by a general 2D scatterer is
$$
P_{abs}=\frac{1}{2}\text{Re}\int_A\mathbf E^*\cdot\mathbf K \text{ dA}=\frac{1}{2}\text{Re}\int_A\mathbf E^*\cdot\mathbf K \text{ dA}
$$


# Scattered Power

The power scattered by a general 2D scatterer is
$$
P_{scat}=\frac{1}{2}\text{Re}\int_A\mathbf E_{inc}^*\cdot\mathbf K \text{ dA}
$$




# Cross Section

Let $\alpha​$ refer to *scattering*, *absorption*, or *extinction*.

## Bounds

According to [Miller et al. 2017](http://pubs.acs.org/doi/10.1021/acs.nanolett.7b02007), assuming a local conductivity relationship between electric field and surface conductivity $\vec K=\bar\sigma\vec E$, the cross section of a 2D scatterer is bound by
$$
\frac{\sigma_\alpha}{A}\leq \beta_\alpha Z_0||\bar\sigma^\dagger(\text{Re}\bar\sigma)^{-1}\bar\sigma||_2
$$
where $\beta_\alpha=1$ for absorption and extinction, and $\beta_\alpha=1/4$ for scattering. Here, $||\cdot||_2​$ is the induced matrix 2-norm (the largest singular value of the matrix).

### Low Frequency Limit

Valid for $\omega\ll \gamma$.
$$
\boxed{\frac{\sigma_{\alpha}}{A}\to\beta_\alpha4\alpha\left(\frac{E_F}{\hbar\gamma}\right),\omega\ll\gamma}
$$
A higher order correction arises from the interband contribution, even below the Fermi level.
$$
\sigma_{\alpha}/A \approx \beta_\alpha4\alpha\left(\frac{E_F}{\hbar\gamma}\right)-\alpha\frac{\hbar\gamma}{E_F}\left[3\left(\omega/\gamma\right)^2-1\right]
$$


### High Frequency Limit

$$
\boxed{\frac{\sigma_\alpha}{A}\to\beta_\alpha\pi\alpha,\omega\gg2E_F/\hbar}
$$



# Polarizibility

See [Garcia de Abajo et al. 2015](https://pubs.rsc.org/en/content/articlelanding/2015/fd/c4fd00216d).

The dipole polarizability of a 2D material sandwiched between two dielectrics $\epsilon_1$ and $\epsilon_2$ is given by $\alpha(\omega)$. It is defined by the induced dipole in response to an applied electric field.
$$
\alpha(\omega)=\frac{1}{E_0}\int\rho^{ind}(\mathbf r) \text{ d}^2\mathbf r
$$
Consider a 2D material with a small characteristic size $D<<\lambda$.
$$
\alpha(\omega)=D^3\sum_{j\in modes} \frac{A_j}{\frac{-1}{n_{eff}^2\eta_j}-\frac{i\omega D}{\sigma(\omega)}}
$$
Here, $n_{eff}=\sqrt{(\epsilon_1+\epsilon_2)/2}​$.

## Sum Rules

### First

For an island of area $A$,
$$
A/D^2 = \sum_j A_j
$$

### Second

$$
\alpha(0)/D^3=-\sum_j \eta_jA_j
$$



# Variables

| Variable                                          | Symbol   | Typical Value                                                | Relations |
| ------------------------------------------------- | -------- | ------------------------------------------------------------ | --------- |
| [Scattering Time](#Scattering Time)               | $\tau$   |                                                              |           |
| [Intrinsic Damping Rate](#Intrinsic Damping Rate) | $\gamma$ |                                                              |           |
| Fermi Level                                       | $E_F$    | $kT$ - 0.6 eV                                                |           |
| Mobility                                          | $\mu$    | 10,000 cm$^2$/V-s [ref](http://www.sciencemag.org/cgi/doi/10.1126/science.1102896) |           |

## Scattering Time

$$
\tau =\mu\hbar\sqrt{n\pi}/ev_F = \mu E_F/ev_F^2
$$

### Examples

The original graphene paper by [Geim and Novoselov](http://www.sciencemag.org/cgi/doi/10.1126/science.1102896) found $\mu=10,000$ cm$^2$/V-s, which yields
$$
\tau =\frac{(1 \text{m}^2/\text{V-s})}{(1\text{eV/V})(1\times 10^6 \text{m/s})^2}E_F=(10^{-12}\text{s/eV})E_F
$$
High quality graphene has $\tau\approx 0.5​$ ps [Woessner](https://www.nature.com/articles/nmat4169).

## Intrinsic Damping Rate

Approximately
$$
\gamma = \tau^{-1}=2\pi^2e^4n_{imp}/\hbar\epsilon_g^2\epsilon
$$
where $n_{imp}$ is the concentration of impurities, $\epsilon_g$ is the permittivity of graphene, and $\epsilon$ is the "typical energy" of electrons. (References 15-17 of Falkovsky et al. 2007)

## Fermi Level

$$
E_F=\hbar v_F\sqrt{\pi n}
$$



## Mobility

$$

$$

