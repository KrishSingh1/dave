# LPJ-GUESS Development Ecosystem

This page provides a high-level overview of how LPJ-GUESS development is organised across version-control systems and teams, and where DAVE now fits.

## Primary development systems

LPJ-GUESS development primarily occurs in **SVN**, on a server hosted by the Lund team.

A large part of the LPJ-GUESS community still works exclusively in SVN, and many contributors do not use Git in day-to-day workflows.

In parallel, the core team maintains a **self-hosted GitLab** server that mirrors trunk and several key branches from the SVN server.

## Typical Git workflow in the LPJ-GUESS community

Among Git users, the common workflow has historically been SVN-like:

- Each developer maintains their own branch on the main shared Git repository.
- Development happens primarily on those per-user branches.
- Collaboration occurs by coordinating changes across those branches.

## DAVE historical workflow

The DAVE team previously followed this same model, with roughly a dozen DAVE branches on the shared upstream Git repository.

## DAVE current workflow

DAVE has now moved to a **fork-based workflow**:

- [A dedicated DAVE fork](https://stormbringer4.nateko.lu.se/dave/lpj-guess-dave) has been created.
- Going forward, all DAVE development is intended to occur on this fork.
- The DAVE documentation website is also maintained in this fork under `docs/`.

## Practical implication

When discussing development workflow for DAVE, this fork should be treated as the canonical development home, while recognising that the broader LPJ-GUESS ecosystem remains anchored in SVN and mirrored into GitLab.
