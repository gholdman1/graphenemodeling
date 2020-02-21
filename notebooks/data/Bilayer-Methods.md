# Bilayer Methods

## Mapping the Dirac Point (Leroy)

[Link to paper](https://aip.scitation.org/doi/10.1063/1.3275755)

Electrostatic equations assume only influence by backgate, and do not assume screening
$$
\text{Dirac Point:  }  E_D=\frac{\hbar^2\pi \epsilon\epsilon_0V_G}{2m^*et_{ox}}
\\
\text{Potential Diff:     }V=\frac{\epsilon d V_G}{2 et_{ox}}
\\
\text{Potential Diff*: } V=-dE=-d(-V_G/t_{ox})...
\\
\text{Dirac Point effective:  }E_D = \frac{\epsilon V_G}{2et_{ox}}\left(\frac{\hbar^2\pi\epsilon_0}{m^*}-d\right)
$$
$V_G$ is the back gate voltage.

Only mention of tip influence

> Also the electric field due to the tip would tend to lessen the effect of the gate and therefore increase the effective mass needed to fit our data.

## Microscopic Polarization (Stroscio)

[Link to paper](http://www.nature.com.ezproxy.library.wisc.edu/articles/nphys1988)

### Extract

They use
$$
E_D = \frac{\hbar^2\pi}{2 m^*}n = \frac{\hbar^2\pi}{2m^*}\alpha(V_G-V_0)
$$
where $\alpha$ comes from parallel plate capacitor equation. They use this to get $m^*$ then and $\gamma_1=m^*v_F^2$ to extract $v_F$.

## Nadj-Perge

[Link to paper](https://export.arxiv.org/ftp/arxiv/papers/1901/1901.02997.pdf)

From supplemental section "Tip-Induced Gating and WF diff":

> In a simple model, the charge density of the TBG underneath the tip can be written as:

$$
n(r,z)=C_{BG}(V_{BG}-V_{Bias})- C_T(r,z)(V_{Bias}-\Delta\Phi)
$$

> ...
>
> By solving the equation for the charge density at charge neutrality for two [current] setpoints [100pA and 1 nA], one gets an estimate for $\Delta\Phi=$ 150-200 mV

This is because $C_T$ changes with a different set point. No screening is mentioned, which makes sense for TBLG.

## Crommie

[Link to paper](http://www.nature.com/articles/nphys1807)

