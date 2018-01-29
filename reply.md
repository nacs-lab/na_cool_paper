# General

Since the last submission, we have made several improvements to our setup and sequence

1. Optimize pulse shaping

    This give us less coupling to unwanted sideband

2. Increase the detuning of the Raman beams from 25GHz to 75GHz

    This allow us to increase the power in the Raman beams and reduce the off-resonance
    scattering rate from them by 2x to 9x.

3. Shorten the cooling sequence

    According to the simulation, some of the pulses we previously use are not
    really necessary and could be causing loss. Removing some pulses groups
    and refining the remainings decrease the number of pulses from 1000 to 540
    and sequence time from 100ms to 53ms.

Combining these improvements, we can now achieve 93.5\% ground state population in 3D with
about 97\% to 98\% per axis which is on par with other Raman sideband cooling result
on neutral atoms. The loss during the cooling sequence is also eliminated completely.

# Referee A

# Referee B

> However, this holds only if the same fidelity for ground state preparation
> as in conventional Raman sideband cooling can be achieved,
> which is not quite the case in the present manuscript.

> The authors report a ground state occupation probability of ~ 90%,
> which is good but a bit below the state-of-the-art for experiments in the Lamb-Dicke regime.
> Could they comment in more detail how much of the shortcoming
> is from experimental limitations and how much is of principal nature?

We believe the limitation are all technical.
Even with the increased detuning, we believe the major limiting factor in the cooling
sequence is still the off-resonance scattering from the Raman beam.
From the numerical simulation,
we can show that without off-resonant scattering from the Raman beams or other technical source
of decoherence, greater than 99% ground state occupation probability is achievable.

> The authors state "The cooling efficiency is limited by spontaneous scattering rate
> (0.5-1 kHz) from the Raman beams, as well as spectral broadening from magnetic field
> fluctuations and trap anharmonicity."
> It would be useful for the reader to detail this quite a bit more.
> In particular in order to enhance the attractiveness to people working
> with different atoms, it would be helpful if the discussed limitations
> could be supported by a analytical or numerical model.

We believe these are the main limitations of our experiement from numerical simulation.
We could not give an universal formula for when these imperfections starts to limit the
cooling performance. However, we have shown that we are not very far from this threshold
and by reducing the Raman beam scattering rate we've seen a big improvement in cooling
performance which agree with what we see in the numerical model.

> How much does the problem of molasses cooling in the optical
> potential (repulsive excited state potential) interfere with the Raman
> cooling.

The effect of the optical potential is eliminated by switching the trap out of phase
with the cooling light. The performance of molasses cooling is therefore not affected
and does not interfere with the Raman cooling.

> The role and origin of the atom loss is not clear to me, however,
> this seems a major roadblock towards achieving large quantum
> registers. What is the origin of the loss and how can it be mitigated?

The loss of the atom was likely caused by the same issues that limits our cooling performance
together with a sequence that is longer than optimal. We have shown that by reducing
the off-resonance scattering and shortening the sequence the atom loss can be eliminated.

# Referee C

> The paper is at most a technical improvement beyond previous works.
> There is no clear reason why the cooling cannot work with large
> vibrational quantum number. A quick survey of some of the previous
> works, I found that Fig. 1 of Ref.[27] already showed how Raman
> cooling can work even in free space.

> Many groups have reported Raman cooling to ground state with Li, K,
> Rb... with even better performance, namely, higher ground state
> fraction. There is no clear reason why sodium is different from the
> others or any experts working on atom tweezers will be excited to
> switch their schemes to the reported one. If such scheme is indeed the
> best for sodium, the paper would only interest sodium atom trappers.
