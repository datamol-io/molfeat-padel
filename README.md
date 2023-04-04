# `MolFeat-Padel`: An extension of Molfeat to support PaDEL-Descriptors

`molfeat-padel` is a plugin to the `molfeat` library that add supports for PaDEL descriptors.

## Usage

The following example show how to use the `molfeat-padel` plugin package automatically when installed. In this example, all three scenarios are valid.

1. initializing the calculator from the plugin package

```python

from molfeat.trans import MoleculeTransformer

from molfeat_padel.calc.padel import PadelDescriptors
trans = MoleculeTransformer(featurizer=PadelDescriptors())
```

2. enable autodiscovery and addition of the `PadelDescriptors` as importable attribute to the entry point group `molfeat.calc`

```python
from molfeat.trans import MoleculeTransformer
from molfeat.plugins import load_registered_plugins
load_registered_plugins(add_submodules=True)

# this is now possible
from molfeat.calc import PadelDescriptors
trans = MoleculeTransformer(featurizer=PadelDescriptors())
```

3. auto discovery of PadelDescriptors 

```python
from molfeat.trans import MoleculeTransformer
import molfeat_padel

trans = MoleculeTransformer(featurizer="PadelDescriptors")
# works because PadelDescriptors is imported in the root init of molfeat_padel
```

```python
from molfeat.trans import MoleculeTransformer
from molfeat.plugins import load_registered_plugins
load_registered_plugins(add_submodules=True)
trans = MoleculeTransformer(featurizer="PadelDescriptors")
```

## Installation

### Conda

Use conda:

```bash
mamba install conda-forge molfeat-padel
```

### Dependencies

The only dependencies of `molfeat-padel` are `padelpy` and `molfeat`

## Changelog
See the latest changelogs at [CHANGELOG.rst](./CHANGELOG.rst).

## Maintainers

- @maclandrol

## License

Under the Apache-2.0 license. See [LICENSE](LICENSE).
