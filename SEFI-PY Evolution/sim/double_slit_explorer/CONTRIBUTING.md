# Contributing

Keep the experience understandable to a curious non-specialist. Place mathematical detail in the advanced view and the physics notes, and keep a plain-language explanation adjacent to each scientific concept.

## Scientific changes

- Identify each statement as established physics, a numerical result under stated assumptions, or a proposed interpretation.
- Support established claims with primary references. Project manuscripts can document proposals but do not establish them experimentally.
- Never present path-bundle animations as tracked particle trajectories, or identity animations as measured observables.
- If a new theory produces probabilities, implement and test its independent prediction separately. Do not silently reuse the quantum reference and label it a new derivation.
- Preserve previous model versions when changing assumptions. Document parameters, units, sampling, and evidence limits.

## Code review

Run `npm test` and `npm run build`. Check every chapter, keyboard controls, small screens, reduced motion, all slider extremes, reset, and image export. Physical parameter changes must clear accumulated samples; display-only changes must not change the detector distribution. Do not add tracking or external AI calls without a documented product decision.

Suggested reviewer questions:

1. Which displayed quantities are calculated, and which are illustrations?
2. Does the path marker recover full visibility at zero distinguishability and remove only the interference term at full distinguishability?
3. Which assumptions are needed for the visibility–distinguishability equality?
4. Does the proposed geometry make any prediction different from the reference model?
5. What mathematical law would be needed to make the DEFI coupling predictive?
6. What observations could falsify a future extension?

These questions may be used with a reviewer's own GitHub Copilot. This app does not expose an author's private Copilot conversations or grant access to their session.
