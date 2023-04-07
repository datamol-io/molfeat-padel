
# `molfeat-padel`: A molfeat plugin that adds support for PaDEL-Descriptors

<div align="center">
    <img src="docs/images/logo-title.svg" width="100%">
</div>

<p align="center">
    <b>molfeat - the hub for all your molecular featurizers </b> <br />
</p>
<p align="center">
  <a href="https://molfeat-docs.datamol.io/" target="_blank">
      Docs
  </a> | 
  <a href="https://molfeat.datamol.io/" target="_blank">
      Homepage
  </a>
</p>

---

[![PyPI](https://img.shields.io/pypi/v/molfeat-padel)](https://pypi.org/project/molfeat-padel/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/molfeat-padel)](https://pypi.org/project/molfeat-padel/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/molfeat-padel)](https://pypi.org/project/molfeat-padel/)
[![license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/datamol-io/molfeat-padel/blob/main/LICENSE)
[![test](https://github.com/datamol-io/molfeat-padel/actions/workflows/test.yml/badge.svg)](https://github.com/datamol-io/molfeat-padel/actions/workflows/test.yml)
[![code-check](https://github.com/datamol-io/molfeat-padel/actions/workflows/code-check.yml/badge.svg)](https://github.com/datamol-io/molfeat-padel/actions/workflows/code-check.yml)
[![release](https://github.com/datamol-io/molfeat-padel/actions/workflows/release.yml/badge.svg)](https://github.com/datamol-io/molfeat-padel/actions/workflows/release.yml)

## Overview

`molfeat-padel` is an extension to `molfeat` that adds support for PaDEL descriptors. 

To learn more about [molfeat](https://github.com/datamol-io/molfeat), please visit https://molfeat.datamol.io/. 

To learn more about the plugin system of molfeat, please see [extending molfeat](https://molfeat-docs.datamol.io/stable/developers/create-plugin.html)

## Installation

You can install `molfeat-padel` with:

```bash
mamba install -c conda-forge molfeat
```
or 

```bash
pip install molfeat-padel
```

`molfeat-padel` depends on `molfeat` and `padelpy`

## Usage

The following example shows how to use the `molfeat-padel` plugin package automatically when installed. All scenarios highlighted in this example are valid:

1. using the `PadelTransformer` directly from `molfeat-padel`

If you intend to use PaDEL descriptors on a list of molecules, this is one of the recommended approach. This approach can significantly improve the efficiency of the process.


```python
from molfeat_padel.trans import PadelTransformer

mol_transf = PadelTransformer()
```

2. interacting directly with `molfeat` by loading the plugins

```python
# put this somewhere in you code
from molfeat.plugins import load_registered_plugins
# in this example we only want to load the molfeat_padel plugin automatically
load_registered_plugins(add_submodules=True, plugins=["molfeat_padel"])
```

```python
# this is now possible
from molfeat.trans import PadelTransformer
mol_transf = PadelTransformer()
```

3. initializing the calculator from the plugin package

```python

from molfeat.trans import MoleculeTransformer

from molfeat_padel.calc import PadelDescriptors
mol_transf = MoleculeTransformer(featurizer=PadelDescriptors())
```

4. enable autodiscovery and addition of the `PadelDescriptors` as importable attribute to the entry point group `molfeat.calc`

```python
# put this somewhere in you code
from molfeat.trans import MoleculeTransformer
from molfeat.plugins import load_registered_plugins
load_registered_plugins(add_submodules=True)
```

```python
# this is now possible
from molfeat.calc import PadelDescriptors
mol_transf = MoleculeTransformer(featurizer=PadelDescriptors())
```

```python
# this is also possible
mol_transf = MoleculeTransformer(featurizer="PadelDescriptors")
```

5. auto discovery of PadelDescriptors 

```python
from molfeat.trans import MoleculeTransformer
import molfeat_padel

mol_transf = MoleculeTransformer(featurizer="PadelDescriptors")
# works because PadelDescriptors is imported in the root init of molfeat_padel
```
### Dependencies

The only dependencies of `molfeat-padel` are `padelpy` and `molfeat`

## Changelog
See the latest changelogs at [CHANGELOG.rst](./CHANGELOG.rst).

## Maintainers

- @maclandrol

## License

Under the Apache-2.0 license. See [LICENSE](LICENSE).
