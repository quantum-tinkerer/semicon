# semicon - k·p simulations made easy

Note: this is work in progress, comments and ideas are more than welcomed!

The goal of this package is to provide easy to use framwerok for performing k·p simulations.
There will be two main components of this framework: parameters and models.

This package is suppose to remove an overhead and a boiler plate that appears when doing the k·p simulations.
It is suppose to make things easier, therefore userfriendly interface is a priority.


# Models

Models are based on the k·p Kane Hamiltonian, symmetrized following Burt-Foreman approach.
Users should be able to choose which components (Zeeman, Dresselhaus, strain, ...) and which bands ('gamma_6c', 'gamma_8v', ...) shall be included.


# Parameters

Parameters bank should include most common semiconducting materials.
Helper functions should handle parameter renormalization based on included bands and way to avoid spurious solutions.


# Helper functions

Goal of this framework is to make k.p simulations easy.
However, building the actual Kwant system and peforming the simulation should be done by user and this library only provide helper tools that are specific to problems encountered in k.p simulations


# Requirements

* This package is suppose to work with latest stable release of Kwant.
* This package will require [SciPy 1.2](https://github.com/scipy/scipy/milestone/36) that is scheduled to be released in November 2018.
This is due to usage of [scipy.spatial.transform.Rotation](https://scipy.github.io/devdocs/generated/scipy.spatial.transform.Rotation.html#scipy.spatial.transform.Rotation) to find the rotation vector of a given rotation matrix.
* Because project is based on discretizer sympy is required, however, due to compatibility [issue](https://gitlab.kwant-project.org/kwant/kwant/issues/225) it must be in version lower than 1.2


# Installation
As this package is pure python standard python installation from sources with
```
python setup.py build
python setup.py install
```
should be working without problems.
The only non-trivial dependency, Kwant, that could cause problemS if not installed beforehand can be either obtained via ``conda`` or any other installation means explained on its [homepage](https://kwant-project.org/).

Direct installation from git is also possible (and favoured as long as conda
package is not out there):
```
pip install git+https://gitlab.kwant-project.org/semicon/semicon.git
```

Note that due to active development master branch may not be stable.
To install version that has been already used in research project use
```
pip install git+https://gitlab.kwant-project.org/semicon/semicon.git@v0.1.0
```


# Tips about developing inside docker container

One can easily use a [rafalskolasinski/science](https://github.com/RafalSkolasinski/science-docker)
Docker container for a development of this project.
Assuming that ``semicon`` folder is ``~/work/semicon`` do:
```
docker pull rafalskolasinski/science:semicon
docker run -d -p 8888:8888 --name semicon \
    -v ~/work/semicon:/src -v ~/work/semicon/notebooks:/home/jovyan/work \
    rafalskolasinski/science:semicon
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
