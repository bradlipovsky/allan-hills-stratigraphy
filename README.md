# Allan Hills cul-de-sac stratigraphy

This project asks what ice-flow and ice-sheet history could place approximately 5 Ma ice directly on top of Pleistocene ice in the Allan Hills cul-de-sac. Older ice also occurs below the Pleistocene unit. The target is therefore a vertically inverted sequence, from top to bottom: old ice, Pleistocene ice, and older ice. It is not merely a lateral boundary between old and young ice. The first objective is not a detailed ice-sheet model. It is to identify the smallest quantitative models that can rule mechanisms in or out.

## Constraints that every hypothesis must satisfy

Let the old and young ice ages be $A_o$ and $A_y$. Let $w_{obs}$ be the measured 10–90% width of the water-isotope step and $w_d$ the physical width inferred after accounting for sample averaging. The observations impose three distinct constraints:

1. The missing-age constraint is $\Delta A=A_o-A_y$, with $A_o\approx5$ Ma and $A_y$ in the Pleistocene. A successful history must remove, bypass, or fold away the intervening ages.
2. The contact-age constraint comes from diffusion. For an initially sharp one-dimensional step with effective diffusivity $D$,

$$
w_d=4\,\text{erf}^{-1}(0.8)\sqrt{Dt}\approx3.62\sqrt{Dt},
\qquad
t_j=\frac{1}{D}\left(\frac{w_d}{3.62}\right)^2.
$$

   Here $t_j$ is the time since the two isotope populations first became adjacent. The existing estimate places $t_j$ at no more than a few tens of thousands of years. Temperature history, strain thinning, sampling resolution, and anisotropic diffusion should be propagated as uncertainties rather than hidden in a single value of $D$.

3. The order-reversal constraint is vertical. For elevation $z_t>z_m>z_b$, a successful section must contain

$$
A(z_t)\approx 5\ \mathrm{Ma},\qquad
A(z_m)=A_y,\qquad
A(z_b)>A_y.
$$

   Thus the Pleistocene unit is a local age minimum bracketed above and below by older ice. A model that only puts old and young packets next to one another horizontally does not explain the observation. The modeled layers must overturn, repeat, or be cut by a thrust or unconformity so that a vertical core recovers the observed order.

These distinctions matter: a mechanism can make a multi-million-year hiatus long ago and still fail because it leaves the two packets in contact for too long, or because it produces lateral juxtaposition without the vertical age reversal.

The time-invariant hypothesis is the negative control. In steady flow, ice age satisfies $\mathbf{u}\cdot\nabla A=1$. A continuous velocity field with ordinary surface and inflow age conditions produces a continuous age field except at an imposed boundary or separatrix. The steady viscous and diagnostic Coulomb calculations in [`ice-wall`](https://github.com/bradlipovsky/ice-wall) did not create the observed vertical order reversal. We will make that no-go result more rigorous only if it becomes necessary; further constitutive complexity is not an initial priority.

This negative control does not exclude passive folding in a steady background flow. [Waddington, Bolzan, and Alley (2001)](https://doi.org/10.3189/172756501781831756) separated folding into two stages: a transient process first displaces a layer from its steady-state orientation, then large-scale flow either overturns that wrinkle by simple shear or flattens it by pure shear. This distinction is relevant to both the Noll hypothesis and the steady-state hypothesis. A steady flow field can amplify a pre-existing wrinkle, but perfectly steady stratigraphy cannot spontaneously reorder itself without a seed disturbance.

For a material layer segment with slope $m=\mathrm{d}z/\mathrm{d}x$, Waddington and colleagues wrote the kinematic evolution as

$$
\frac{\mathrm{D}m}{\mathrm{D}t}
=-m^2\frac{\partial u}{\partial z}
-m\left(\frac{\partial u}{\partial x}-\frac{\partial w}{\partial z}\right)
+\frac{\partial w}{\partial x}.
$$

The first term contains the vertical shear that can overturn a sloping limb; the second represents the pure-shear tendency to flatten it. Their dimensionless shear number,

$$
S=\frac{\partial u/\partial z+\partial w/\partial x}
        {\partial u/\partial x-\partial w/\partial z},
$$

compares these effects. Relative to a local steady isochrone of slope $m_0$, the approximate critical disturbance is $\Delta m_{crit}\approx-S^{-1}-2m_0$. Because strain rates vary along a path, the stronger test is to integrate the deformation gradient and ask whether a seeded segment actually passes through vertical and overturns. Every model below should therefore report both how a wrinkle could be seeded and whether its subsequent finite strain creates an old–Pleistocene–older vertical section.

## Possible histories

### 1. Recent collision with an old cul-de-sac reservoir

Ancient ice remains nearly stagnant and cold-based in a topographic pocket. A change in ice thickness, flow direction, divide position, or upstream flux later routes a lobe of Pleistocene ice against or beneath that reservoir. The required outcome is a thrust-like or recumbent geometry in which the old reservoir lies above the incoming Pleistocene ice while older ice remains below. The structural analogue is a surge front running into stagnant ice. At Variegated Glacier, this geometry produced intense shortening, folding, and faulting ([Raymond et al., 1987](https://doi.org/10.1029/JB092iB09p09037)); thrusts dipping about 40° up-glacier accounted for as much as 50% of local shortening ([Moore et al., 2010](https://doi.org/10.1029/2009JF001307)). Allan Hills is not a temperate surge glacier, so this analogy concerns the kinematics of a moving compressional front, not its basal hydrology.

The minimal model is a two-dimensional full-Stokes flow band with prescribed transient boundary forcing. Track age as a material label using

$$
\frac{\partial A}{\partial t}+\mathbf{u}\cdot\nabla A=1,
$$

and track the isotope step without numerical diffusion. Vary the event age before present $t_e$, duration $T_e$, speed-up $S=U_e/U_0$, front speed $c_f$, integrated shortening $\int-\dot\epsilon_{xx}\,dt$, flow-direction change, and initial separation $L_{sep}$. A necessary kinematic condition is $U_eT_e\gtrsim L_{sep}$, while the final contact must satisfy $t_e\lesssim t_j$. The model should test whether old ice can be carried over the young packet and leave a recumbent fold, repeated section, thrust, or shear zone that yields old–Pleistocene–older in a vertical core. This is the highest-priority world because it requires time dependence but not an exotic rheology.

### 2. Surface-scour unconformity followed by recent deformation

Persistent wind scour can create an unconformity by preventing accumulation while ice passes through a fixed surface zone. Radar observations south of Dome A trace such unconformities for hundreds of kilometers and through much of the ice column; kinematic models are consistent with stable scour and flow over multiple glacial cycles ([Fudge et al., 2026](https://doi.org/10.1017/jog.2026.10171)).

A scour-only history is viable only if the time since deposition of the adjoining young layer also satisfies the diffusion limit. Otherwise the useful hypothesis is composite: scour first creates a missing-age boundary, and a later fold, truncation, or flow reorganization places the old unit above the young one. Model the scour zone as $a(x,t)=0$ or $a<0$, then test whether its downstream unconformity can be folded or repeated into an old–Pleistocene–older vertical section without diffusing the isotope step for too long.

### 3. Traveling basal slippery patch

A moving basal traction anomaly can produce uplift, subsidence, overturned layers, and thickness-scale folds. In the thermomechanical experiments of [Wolovick et al. (2014)](https://doi.org/10.1002/2014GL062248), a small water-flux perturbation generated traveling slippery patches, uplifted basal ice to about half the ice thickness, and produced 20–50% flux increases lasting roughly 1 kyr with onsets near 100 yr.

The first test should prescribe a moving low-drag patch rather than reproduce the full thermal feedback. Its parameters are patch width $L_p/H$, speed $c_p/U$, traction contrast, and lifetime. Success requires an overturned or repeated age sequence that places old ice above Pleistocene ice while retaining older ice below. The Waddington criterion provides a direct diagnostic of whether the patch merely tilts layers or seeds limbs that background shear can overturn. Only then should we test whether geothermal, frictional, and conductive heat budgets permit the required basal state at Allan Hills. A warm-based mechanism that cannot pass that thermal check is rejected even if its kinematics work.

### 4. Migrating Ice Wall

The Ice Wall is a steep surface feature that may be a leeward snow drift accumulated over approximately 100 kyr. The steady calculations in [`ice-wall`](https://github.com/bradlipovsky/ice-wall) tested flow beneath a fixed wall. If the accumulation pattern and wall migrate, however, the resulting traveling surface-slope anomaly could generate a moving zone of convergence and vertical motion. It is the surface analogue of a traveling basal slippery patch and may be coupled to surface scour: erosion on one side and deposition on the other would translate both the wall and an unconformity.

Compare a fixed wall with a migrating surface mass-balance dipole or prescribed translating surface feature. Vary wall height $\Delta H/H$, width $L_w/H$, migration speed $c_w/U$, build time $T_w$, and accumulation–ablation amplitude. Track particles and isochrones to determine whether the moving wall can seed a fold or truncate the age field so that old ice overlies Pleistocene ice within $t_j$, with older ice beneath both. Diagnose the layer slopes and finite strain using the Waddington criterion. A successful model must also predict a wall-migration distance $c_wT_w$, layer geometry, and surface mass-balance pattern that can be checked with radar and field observations.

These histories are not mutually exclusive. The likely composite world is long-term preservation of old ice, formation or advection of an age unconformity or wrinkle, and a much more recent event that creates the present sharp inverted contact. In particular, a migrating Ice Wall, surface scour, or a basal patch could seed a disturbance that later steady shear overturns, as envisaged by the Waddington framework.

## Initial two-dimensional full-Stokes tests

Issue [#11](https://github.com/bradlipovsky/allan-hills-stratigraphy/issues/11) replaces the reduced collision, basal-patch, and Ice-Wall calculations with deliberately simple full-Stokes tests. Each notebook solves for horizontal velocity, vertical velocity, and pressure using a P2–P1 mixed finite-element pair, the complete symmetric velocity gradient, incompressibility, and gravity. Initially ordered material layers are then advected through the solved velocity and sampled with objective vertical synthetic cores.

| Scenario | Idealized forcing | Maximum resolved vertical shear | Result |
|---|---|---:|---|
| [Collision](notebooks/full_stokes_collision.ipynb) | Opposing side inflow for 30 yr; compensating stress-free surface outflow | $6.00\times10^{-2}$ yr$^{-1}$ | 0 overturned segments; no bracketed Pleistocene minimum |
| [Traveling basal patch](notebooks/full_stokes_basal_patch.ipynb) | 8 km low-drag Robin patch; 600 yr in a frame moving at 40 m yr$^{-1}$ | $2.98\times10^{-3}$ yr$^{-1}$ | 0 overturned segments; no bracketed Pleistocene minimum |
| [Migrating Ice Wall](notebooks/full_stokes_ice_wall.ipynb) | 200 m surface step; 1 km of relative migration over 100 kyr | $6.62\times10^{-6}$ yr$^{-1}$ in the wall zone | 0 overturned segments; no bracketed Pleistocene minimum |

All three conclusions persist when both mesh dimensions are doubled, and their boundary mass-flux residuals are far below 1%. These are negative results for the stated smooth, linear-viscous end members—not proofs that the broader histories are impossible. In particular, the collision calculation requires a large compensating surface flux, the basal patch is prescribed rather than thermomechanically generated, and the Ice Wall is a fixed geometry viewed in a translating frame. The calculations nevertheless show that resolved depth-dependent shear alone does not produce the observed old–Pleistocene–older order: the forcing must first create a sufficiently steep material wrinkle, thrust, or unconformity for subsequent shear to overturn.

## Minimal modeling sequence and success metrics

1. Assemble the measured ages, isotope profiles, sample-response length, core locations, contact orientation, ice thickness, bed and surface geometry, accumulation/ablation, and present velocity. Record uncertainties and do not interpolate across the contact.
2. Recalculate the diffusion bound by fitting an error-function step convolved with the sampling kernel. Report distributions for $w_d$ and $t_j$, not only best values.
3. Use simple two-dimensional full-Stokes models for the collision, basal-patch, and Ice-Wall worlds, while treating surface scour as a boundary-history problem. Prescribe the transient forcing when possible, and keep age and isotope tracers Lagrangian or demonstrate convergence so that numerical diffusion is smaller than the observed width.
4. Add thermomechanical feedback or basal hydrology only after an idealized mechanical model passes the observational tests.

A scenario passes the initial investigation if it satisfies all of the following:

- Its modeled old–young age difference overlaps the dating intervals for $A_o$ and $A_y$.
- A vertical synthetic core contains a Pleistocene age minimum bracketed by older ice: old–Pleistocene–older from top to bottom. Report the depths, thicknesses, and both contact orientations. Lateral contact alone is a failure.
- After physical diffusion and the same measurement convolution, its modeled 10–90% isotope width overlaps $w_{obs}$; equivalently, direct juxtaposition occurs no earlier than the inferred $t_j$.
- Contact position and orientation fit the available core or radar observations within one stated observational resolution, or have normalized RMS misfit no greater than one when formal uncertainties exist.
- Ice mass is conserved to 1% in the toy model, and the required displacement, thickness change, velocity, accumulation, basal temperature, and heat budget remain within observationally defensible ranges.
- It predicts at least one independent sign test: continuation of the contact, a repeated or overturned sequence, a fold or shear zone, or wall-related layer geometry.
- The conclusion is stable to a factor-of-two resolution change and to a documented parameter sweep over the uncertain inputs.

The project succeeds when these tests rank the possible histories and identify the most discriminating next observation. A unique reconstruction of the Allan Hills geometry is not required.

## Starting literature

- [Shackleton et al. (2025), Miocene and Pliocene ice and air from the Allan Hills blue ice area](https://doi.org/10.1073/pnas.2502681122)
- [Higgins et al. (2015), Atmospheric composition 1 million years ago from blue ice in the Allan Hills](https://doi.org/10.1073/pnas.1420232112)
- [Kehrl et al. (2018), Evaluating the duration and continuity of potential climate records from the Allan Hills](https://doi.org/10.1029/2018GL077511)
- [Fudge et al. (2026), Stability of interior East Antarctic wind scour and ice flow for multiple glacial-interglacial cycles](https://doi.org/10.1017/jog.2026.10171)
- [Wolovick et al. (2014), Traveling slippery patches produce thickness-scale folds in ice sheets](https://doi.org/10.1002/2014GL062248)
- [Waddington et al. (2001), Potential for stratigraphic folding near ice-sheet centers](https://doi.org/10.3189/172756501781831756)
- [Jacobson and Waddington (2004), Recumbent folding in ice sheets: a core-referential study](https://doi.org/10.3189/172756504781830204)
- [Raymond (1987), How do glaciers surge?](https://doi.org/10.1029/JB092iB09p09121)
- [Raymond et al. (1987), Propagation of a glacier surge into stagnant ice](https://doi.org/10.1029/JB092iB09p09037)
