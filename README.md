# Allan Hills cul-de-sac stratigraphy

This project asks what ice-flow and ice-sheet history could place approximately 5 Ma ice directly beside Pleistocene ice in the Allan Hills cul-de-sac. The first objective is not a detailed ice-sheet model. It is to identify the smallest quantitative models that can rule mechanisms in or out.

## Constraints that every hypothesis must satisfy

Let the old and young ice ages be $A_o$ and $A_y$. Let $w_{obs}$ be the measured 10–90% width of the water-isotope step and $w_d$ the physical width inferred after accounting for sample averaging. The observations impose two distinct constraints:

1. The missing-age constraint is $\Delta A=A_o-A_y$, with $A_o\approx5$ Ma and $A_y$ in the Pleistocene. A successful history must remove, bypass, or fold away the intervening ages.
2. The contact-age constraint comes from diffusion. For an initially sharp one-dimensional step with effective diffusivity $D$,

$$
w_d=4\,\text{erf}^{-1}(0.8)\sqrt{Dt}\approx3.62\sqrt{Dt},
\qquad
t_j=\frac{1}{D}\left(\frac{w_d}{3.62}\right)^2.
$$

   Here $t_j$ is the time since the two isotope populations first became adjacent. The existing estimate places $t_j$ at no more than a few tens of thousands of years. Temperature history, strain thinning, sampling resolution, and anisotropic diffusion should be propagated as uncertainties rather than hidden in a single value of $D$.

The distinction matters: a mechanism can make a multi-million-year hiatus long ago and still fail because it leaves the two packets in contact for too long.

The time-invariant hypothesis is the negative control. In steady flow, ice age satisfies $\mathbf{u}\cdot\nabla A=1$. A continuous velocity field with ordinary surface and inflow age conditions produces a continuous age field except at an imposed boundary or separatrix. The steady viscous and diagnostic Coulomb calculations in [`ice-wall`](https://github.com/bradlipovsky/ice-wall) did not create the observed juxtaposition. We will make that no-go result more rigorous only if it becomes necessary; further constitutive complexity is not an initial priority.

## Possible histories

### 1. Recent collision with an old cul-de-sac reservoir

Ancient ice remains nearly stagnant and cold-based in a topographic pocket. A change in ice thickness, flow direction, divide position, or upstream flux later routes Pleistocene ice against or over that reservoir. This directly separates the preservation problem from the recent-contact problem.

The minimal model is a two-dimensional incompressible flow band with prescribed transient velocity. Track age as a material label using

$$
\frac{\partial A}{\partial t}+\mathbf{u}\cdot\nabla A=1,
$$

and track the isotope step without numerical diffusion. Vary the event age before present $t_e$, duration $T_e$, speed-up $S=U_e/U_0$, flow-direction change, and initial separation $L_{sep}$. A necessary kinematic condition is $U_eT_e\gtrsim L_{sep}$, while the final contact must satisfy $t_e\lesssim t_j$. This is the highest-priority world because it requires time dependence but not an exotic rheology.

### 2. Surface-scour unconformity followed by recent deformation

Persistent wind scour can create an unconformity by preventing accumulation while ice passes through a fixed surface zone. Radar observations south of Dome A trace such unconformities for hundreds of kilometers and through much of the ice column; kinematic models are consistent with stable scour and flow over multiple glacial cycles ([Fudge et al., 2026](https://doi.org/10.1017/jog.2026.10171)).

A scour-only history is viable only if the time since deposition of the adjoining young layer also satisfies the diffusion limit. Otherwise the useful hypothesis is composite: scour first creates a missing-age boundary, and a later fold, truncation, or flow reorganization brings the measured packets into contact. Model the scour zone as $a(x,t)=0$ or $a<0$, then test whether its downstream unconformity can be folded into the observed orientation without diffusing the isotope step for too long.

### 3. Traveling basal slippery patch

A moving basal traction anomaly can produce uplift, subsidence, overturned layers, and thickness-scale folds. In the thermomechanical experiments of [Wolovick et al. (2014)](https://doi.org/10.1002/2014GL062248), a small water-flux perturbation generated traveling slippery patches, uplifted basal ice to about half the ice thickness, and produced 20–50% flux increases lasting roughly 1 kyr with onsets near 100 yr.

The first test should prescribe a moving low-drag patch rather than reproduce the full thermal feedback. Its parameters are patch width $L_p/H$, speed $c_p/U$, traction contrast, and lifetime. Success requires an overturned or repeated age sequence with the correct contact geometry. Only then should we test whether geothermal, frictional, and conductive heat budgets permit the required basal state at Allan Hills. A warm-based mechanism that cannot pass that thermal check is rejected even if its kinematics work.

### 4. Surge-front compression into stagnant ice

The structural analogue is a short-lived fast-flow lobe running into a stagnant toe. At Variegated Glacier, a surge propagated into stagnant ice and produced intense shortening, folding, and faulting ([Raymond et al., 1987](https://doi.org/10.1029/JB092iB09p09037)). Thrusts dipping about 40° up-glacier accounted for as much as 50% of local shortening in the observed temperate-glacier example ([Moore et al., 2010](https://doi.org/10.1029/2009JF001307)).

Allan Hills is not a temperate surge glacier, so basal hydrology should not be imported as an assumption. Test the geometry first with a moving compressional front described by front speed $c_f$, speed ratio $S$, duration $T_e$, and integrated shortening $\int-\dot\epsilon_{xx}\,dt$. The hypothesis is supported only if a plausible transient can vault young ice onto an old stagnant packet and leave a fold, repeated section, or shear zone that radar, fabric, or borehole data could detect.

### 5. Basal accretion or selective removal

Basal freeze-on can raise old or debris-bearing ice, while basal melting or surface ablation can selectively remove intermediate ages. This world predicts provenance changes—gas-poor ice, debris, unusual chemistry or crystallographic fabric—rather than merely a clean meteoric-ice contact. It can be screened with existing core observations before modeling. If both packets are ordinary meteoric ice and no removal surface is present, this mechanism becomes a component of another history rather than the primary explanation.

These histories are not mutually exclusive. The likely composite world is long-term preservation of old ice, formation or advection of an age unconformity, and a much more recent event that creates the present sharp contact.

## Minimal modeling sequence and success metrics

1. Assemble the measured ages, isotope profiles, sample-response length, core locations, contact orientation, ice thickness, bed and surface geometry, accumulation/ablation, and present velocity. Record uncertainties and do not interpolate across the contact.
2. Recalculate the diffusion bound by fitting an error-function step convolved with the sampling kernel. Report distributions for $w_d$ and $t_j$, not only best values.
3. Use prescribed, mass-conserving two-dimensional velocity histories to test worlds 1–4. Keep age and isotope tracers Lagrangian or demonstrate convergence so that numerical diffusion is smaller than the observed width.
4. Add stress, temperature, or basal hydrology only after a kinematic history passes the observational tests.

A scenario passes the initial investigation if it satisfies all of the following:

- Its modeled old–young age difference overlaps the dating intervals for $A_o$ and $A_y$.
- After physical diffusion and the same measurement convolution, its modeled 10–90% isotope width overlaps $w_{obs}$; equivalently, direct juxtaposition occurs no earlier than the inferred $t_j$.
- Contact position and orientation fit the available core or radar observations within one stated observational resolution, or have normalized RMS misfit no greater than one when formal uncertainties exist.
- Ice mass is conserved to 1% in the toy model, and the required displacement, thickness change, velocity, accumulation, basal temperature, and heat budget remain within observationally defensible ranges.
- It predicts at least one independent sign test: continuation of the contact, a repeated or overturned sequence, a fold or shear zone, or a provenance/fabric change.
- The conclusion is stable to a factor-of-two resolution change and to a documented parameter sweep over the uncertain inputs.

The project succeeds when these tests rank the possible histories and identify the most discriminating next observation. A unique full-Stokes reconstruction is not required.

## Starting literature

- [Shackleton et al. (2025), Miocene and Pliocene ice and air from the Allan Hills blue ice area](https://doi.org/10.1073/pnas.2502681122)
- [Higgins et al. (2015), Atmospheric composition 1 million years ago from blue ice in the Allan Hills](https://doi.org/10.1073/pnas.1420232112)
- [Kehrl et al. (2018), Evaluating the duration and continuity of potential climate records from the Allan Hills](https://doi.org/10.1029/2018GL077511)
- [Fudge et al. (2026), Stability of interior East Antarctic wind scour and ice flow for multiple glacial-interglacial cycles](https://doi.org/10.1017/jog.2026.10171)
- [Wolovick et al. (2014), Traveling slippery patches produce thickness-scale folds in ice sheets](https://doi.org/10.1002/2014GL062248)
- [Raymond (1987), How do glaciers surge?](https://doi.org/10.1029/JB092iB09p09121)
- [Raymond et al. (1987), Propagation of a glacier surge into stagnant ice](https://doi.org/10.1029/JB092iB09p09037)
- [Bell et al. (2011), Widespread persistent thickening of the East Antarctic Ice Sheet by freezing from the base](https://doi.org/10.1126/science.1200109)
