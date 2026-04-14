---
name: settings-py
description: "Use when: configuring HPC clusters in settings.py, adding new clusters, modifying Q/Schrodinger paths, adjusting forcefield registry. Always check settings.py after qligfep_init.sh."
applyTo: "**/settings.py"
---

# settings.py Configuration Guide

## Overview

`settings.py` is populated by `qligfep_init.sh` during initial setup. It contains:
- **Path configurations**: Q software, Schrodinger, forcefields
- **HPC cluster profiles**: Submit job scripts, module loads, resource specs
- **Forcefield registry**: Available FFs for converters and CLI validation
- **Template locations**: Q input files for MD/FEP/equilibration

---

## Key Sections

### 1. Path Configuration (Set by `qligfep_init.sh`)

```python
# Q Software installation directory
Q_DIR = '/path/to/Q6'

# Schrodinger installation (for PyMOL, ffld_server, Protein Prep)
SCHROD = '/path/to/schrodinger'

# Directories (relative paths from repo root)
FF_DIR = os.path.join(QLIGFEP, 'FF')
INPUT_DIR = os.path.join(QLIGFEP, 'INPUTS')
```

**Edit when**:
- Q software installed in new location
- Schrodinger moved or updated
- Template files reorganized

---

### 2. Forcefield Registry

```python
FORCEFIELDS = {
    'OPLS2015': 'OPLS2015.prm',
    'OPLS2005': 'OPLS2005.prm',
    'AMBER14sb': 'AMBER14sb.prm',
    'CHARMM36': 'CHARMM36.prm',
}
```

**Edit when**:
- Adding new FF variant (see `/add-forcefield` prompt)
- Renaming FF file in `FF/` directory
- Deprecating old FF versions

**Format**: `'FF_NAME': 'filename.prm'`

---

### 3. HPC Cluster Configuration (Most Common Edit)

Each cluster is a dict with execution commands and resource specs:

```python
CLUSTERS = {
    'localhost': {
        'Q_DIR': Q_DIR,
        'QDYN': os.path.join(Q_DIR, 'bin', 'qdyn'),
        'QPREP': os.path.join(Q_DIR, 'bin', 'qprep'),
        'QFEP': os.path.join(Q_DIR, 'bin', 'qfep'),
        'QCALC': os.path.join(Q_DIR, 'bin', 'qcalc'),
        'SUBMIT_CMD': 'bash',  # Local execution
        'MODULES': '',  # No module loads
        'NODES': 1,
        'NTASKS': 1,
        'CPUS_PER_TASK': 1,
        'TIME': '10:00:00',  # HH:MM:SS
        'PARTITION': 'normal',
    },
    'your_cluster': {
        'Q_DIR': '/path/to/Q6/on/cluster',
        'QDYN': '/path/to/Q6/bin/qdyn',
        'QPREP': '/path/to/Q6/bin/qprep',
        'QFEP': '/path/to/Q6/bin/qfep',
        'QCALC': '/path/to/Q6/bin/qcalc',
        'SUBMIT_CMD': 'sbatch',  # SLURM submission
        'MODULES': 'module load intel/2021.3 Q/6.1',
        'NODES': 1,
        'NTASKS': 8,
        'CPUS_PER_TASK': 2,
        'TIME': '24:00:00',
        'PARTITION': 'gpu',  # If GPU cluster
    },
}
```

### Adding a New Cluster

1. **Determine submission system** (SLURM, PBS, SGE, local):
   ```python
   # SLURM: 'sbatch'
   # PBS:   'qsub'
   # SGE:   'qsub'
   # Local: 'bash'
   'SUBMIT_CMD': 'sbatch',
   ```

2. **Get module load commands** from cluster documentation:
   ```bash
   # Example: ssh cluster_login@hostname
   # module avail  # List available modules
   # module list   # Current modules
   ```

3. **Query resource specs**:
   ```bash
   # Get default node/task counts from cluster docs or:
   sinfo -n node1 -o "%N %c %m"  # CPUs, memory
   ```

4. **Test Q installation path** on cluster:
   ```bash
   ssh cluster_login@hostname
   which qdyn  # or: ls /path/to/Q6/bin/
   ```

5. **Add to CLUSTERS dict**:
   ```python
   CLUSTERS = {
       ...
       'my_cluster': {
           'Q_DIR': '/cluster/Q6',
           'QDYN': '/cluster/Q6/bin/qdyn',
           'QPREP': '/cluster/Q6/bin/qprep',
           'QFEP': '/cluster/Q6/bin/qfep',
           'QCALC': '/cluster/Q6/bin/qcalc',
           'SUBMIT_CMD': 'sbatch',
           'MODULES': 'module load intel/2021 openmpi/4.0',
           'NODES': 1,
           'NTASKS': 16,
           'CPUS_PER_TASK': 1,  # Or higher if hyperthreading
           'TIME': '48:00:00',
           'PARTITION': 'normal',
       },
   }
   ```

---

### 4. Default Cluster (Set by `qligfep_init.sh`)

```python
DEFAULT = 'localhost'  # Fallback if -C not specified in CLI
```

**Usage**:
```bash
python QligFEP.py -l1 lig1.pdb -l2 lig2.pdb -f OPLS2015 -S prot.pdb
# Uses DEFAULT cluster from settings.py

python QligFEP.py -l1 lig1.pdb -l2 lig2.pdb -f OPLS2015 -S prot.pdb -C my_cluster
# Uses 'my_cluster' from CLUSTERS dict
```

---

## Validation & Troubleshooting

### Check Runtime Configuration

```python
# Add to settings.py to debug:
import pprint
print("Active Clusters:")
pprint.pprint(CLUSTERS)
print(f"\nDefault: {DEFAULT}")
print(f"Q_DIR: {Q_DIR}")
print(f"FF_DIR: {FF_DIR}")
```

### Common Issues

| Issue | Solution |
|-------|----------|
| "Cluster not found" | Check cluster name in CLUSTERS dict matches `-C` flag |
| Q commands not found | Verify Q_DIR and individual cmd paths in cluster config |
| Module load fails | Check MODULES string syntax for your scheduler (SLURM vs PBS vs SGE) |
| Wrong Q executable | Ensure QDYN/QPREP/QFEP point to correct binary paths on cluster |
| Jobs timeout | Increase TIME limit (HH:MM:SS format) in cluster config |
| Resource exhaustion | Increase NTASKS or CPUS_PER_TASK; check cluster max per user |

### Manual Verification on Cluster

```bash
# SSH to cluster
ssh login@cluster_host

# Load modules manually
module load intel/2021 Q/6.1

# Test each executable
qdyn
qprep
qfep
qcalc

# Verify paths
which qdyn  # Should match QDYN path in settings.py
```

---

## Template Files Reference

Q input templates in `INPUTS/` are formatted with parameters and saved per-system:

```python
# These are auto-loaded from INPUT_DIR:
TEMPLATES = {
    'eq1': 'eq1.inp',
    'eq2': 'eq2.inp',
    'eq3': 'eq3.inp',
    'eq4': 'eq4.inp',
    'eq5': 'eq5.inp',
    'md_0000_1000': 'md_0000_1000.inp',
    'qprep': 'qprep_QligFEP.inp',
    'qfep': 'qfep.inp',
}
```

---

## Post-Init Customization Workflow

1. Run initialization: `bash qligfep_init.sh`
2. Test with local cluster: `python QligFEP.py ... -C localhost`
3. Add new clusters to `CLUSTERS` dict
4. Test each cluster: `python QligFEP.py ... -C cluster_name`
5. Update DEFAULT if needed: `DEFAULT = 'cluster_name'`

