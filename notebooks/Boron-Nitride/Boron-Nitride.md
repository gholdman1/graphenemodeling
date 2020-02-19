[TOC]

# Phonons

Hexagonal boron nitride has an *in-plane* phonon at 1370 cm$^{-1}$ and an out-of-plane phonon at 783 cm$^{-1}$both in the IR. [Geick 1966](https://link.aps.org/doi/10.1103/PhysRev.146.543)

# Permittivity



For notes on modeling in optical simulations, see [Modeling/Optical](#Optical)

# Modeling

## Optical

From [S1 of Brar 2014](https://doi.org/10.1021/nl501096s):

Modeled in a Lorentz oscillator model as
$$
\epsilon_{BN}(\omega)=\epsilon_{\infty}+\sum_i\frac{s_i^2}{\omega_i^2-\omega^2+i\omega\gamma_i}
$$

| Parameter         | Value                | Description                               |
| ----------------- | -------------------- | ----------------------------------------- |
| $\epsilon_\infty$ | 4.95                 | Optical Dielectric Constant (Ref 4 in S1) |
| $\omega_1$        | 1370 cm$^{-1}$       | In-plane phonon frequency                 |
| $s_1^2$           | 3.9x10$^6$ cm$^{-2}$ |                                           |
| $\gamma_1$        | 19 cm$^{-1}$         | Loss                                      |

