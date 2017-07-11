# semicon - a k.p simulations made easy

Note: this is work in progress, comments and ideas are more than welcomed!

The goal of this package is to provide useful tools for doing k.p simulations.
This tools consist of model definitions, material parameters, and helper
functions to generate Kwant template systems.

This package is suppose to remove an overhead and a boiler plate that appears
when doing the k.p simulations. It is suppose to make things easier, therefore
userfriendly interface is a priority.


# Models

Models are based on the k.p Kane Hamiltonian, symmetrized following Burt-Foreman
approach. User should be able to choose which components he wants to include
(Zeeman term, strain, Dresselhaus SOI, etc.) and which bands (standard 8x8 and 14x14 versions).


# Parameters

Parameter bank should contain all standard semiconductor materials, together
with references. Function to renormalize parameters, depending on bands included
in the model, should be also included.


# Helper functions

Helper functions should generate Kwant template system or provide all necessary
information to do so. Building the system and making the simulation itself
should be done by user as it should be a standard Kwant operation.

It is possible that some functions to aid further simulation will be included,
if they will not suit into ``kwant.continuum`` module.



# Requirements

Package (stable release / master branch) will be meant to work with latest
stable release of Kwant. Features that will depend on unreleased functionality
(kwant's master branch) will be provided in the development (unstable) branch
or release.

Because project is based on discretizer sympy is required.
Additional requirments are possible.


# Tips about developing inside docker container

One can use ``kwant-devenv`` for a development of this project.
Assuming that ``semicon`` folder is ``~/projects/semicon`` do:
```
docker pull rafalskolasinski/kwant-devenv
docker run -d -p 8888:8888 --name semicon \
    -v ~/projects/semicon:/src -v ~/projects/semicon/notebooks:/notebooks \
    rafalskolasinski/kwant-devenv
```

This will mount source code in ``/src`` and project notebooks in ``/notebooks``
inside the containier. It will also start ``jupyter notebook`` server running
at ``localhost:8888``.

You can now use ``docker exec semicon build`` and ``docker exec semicon test``
to build the package or run the tests respectively.

You can enter bash inside the container by running
```
docker exec -it semicon bash
```
