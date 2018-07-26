# semicon - k·p simulations made easy

Note: this is work in progress, comments and ideas are more than welcomed!

The goal of this package is to provide useful tools for doing k·p simulations.
This tools consist of model definitions, material parameters, and helper
functions to generate Kwant template systems.

This package is suppose to remove an overhead and a boiler plate that appears
when doing the k·p simulations. It is suppose to make things easier, therefore
userfriendly interface is a priority.


# Models

Models are based on the k·p Kane Hamiltonian, symmetrized following Burt-Foreman
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

# Installation
As this package is pure python standard python installation from sources with
```
python setup.py build
python setup.py install
```
should be working without problems. The only non-trivial dependency, Kwant, that
could cause problem if not installed beforehand can be either obtained via
``conda`` or any other installation means explained on
its [homepage](https://kwant-project.org/).

Direct installation from git is also possible (and favoured as long as conda
package is not out there):
```
pip install git+https://gitlab.kwant-project.org/r-j-skolasinski/semicon.git
```


# Tips about developing inside docker container

One can easily use a [rafalskolasinski/science](https://github.com/RafalSkolasinski/science-docker) 
Docker container for a development of this project.
Assuming that ``semicon`` folder is ``~/projects/semicon`` do:
```
docker pull rafalskolasinski/science
docker run -d -p 8888:8888 --name semicon \
    -v ~/work/semicon:/src -v ~/work/semicon/notebooks:/home/jovyan/work \
    rafalskolasinski/science
```

This will mount source code in ``/src`` and project notebooks in ``~/work``
inside the containier. It will also start ``jupyter notebook`` server running
at ``localhost:8888``. You will need to read jupyter's server token with 
``docker logs semicon`` to access the server.

You can now use ``docker exec semicon build`` and ``docker exec semicon test``
to build the package or run the tests respectively.

You can enter bash inside the container by running
```
docker exec -it semicon bash
```

Nicely formatted output of tests (colors):
```
docker exec -it semicon test -v
```
