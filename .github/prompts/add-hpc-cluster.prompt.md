---
name: add-hpc-cluster
description: "Use when: adding a new HPC cluster to settings.py, configuring SLURM/PBS/SGE clusters, updating Q paths, troubleshooting job submission failures."
---

# Add an HPC Cluster to settings.py

Follow this workflow to integrate a new compute cluster for QligFEP/QresFEP job submission.

## Step 1: Gather Cluster Information

### Connect to cluster and inspect system

```bash
# SSH to cluster login node
ssh username@cluster_host

# Check scheduler system
which sbatch      # SLURM
qstat -a          # PBS
qhost             # SGE

# Verify Q software installation
which qdyn qprep qfep qcalc
# Or check explicit path
ls /cluster/software/Q6/bin/

# Get available modules (if using module system)
module avail
module list
```

### Determine cluster resource specs

```bash
# SLURM: Check node/task config
sinfo -n node1 -o "%N %c %m"     # CPUs, memory per node

# PBS: Check scheduler config
pbsnode -a node1                  # Node resources

# Check your user limits
scontrol show user=yourname       # SLURM
qmgr -c "list server"             # PBS
```

### Find scratch/work directories

```bash
# Most clusters have local scratch for fast I/O
echo $TMPDIR          # Local temp directory
ls -lh $SCRATCH       # Scratch directory path
# Note: store large I/O files here, not $HOME
```

---

## Step 2: Determine Submission System

QligFEP uses cluster-specific submission commands. Identify your cluster system:

| Cluster Type | Submit Command | Module System | Example Location |
|--------------|---|---|---|
| **SLURM** | `sbatch` | Usually: `module load intel openmpi` | High Performance Computing (HPC) clusters |
| **PBS/Torque** | `qsub` | Usually: `module load intel mvapich2` | Older university HPC |
| **Grid Engine (SGE)** | `qsub` | Usually: `module load sge` | Research computing |
| **Local machine** | `bash` | None | Laptop/workstation (debugging) |

---

## Step 3: Configure Module Loads

Most clusters require loading specific software before running jobs.

```bash
# Test module loads interactively
module purge                          # Start fresh
module load intel/2021.3              # Compiler
module load openmpi/4.0.1             # MPI
module list                           # Verify

# Test Q software after modules
qdyn -h                               # Should work
```

**Save the exact module names** for your `settings.py` configuration:
```python
'MODULES': 'module load intel/2021.3 openmpi/4.0.1',
```

---

## Step 4: Test Q Installation Path

```bash
# After loading modules, find exact Q paths
which qdyn
which qprep
which qfep
which qcalc

# Or check a known installation dir
ls -la /cluster/software/Q6/bin/
```

**Common Q installation locations**:
- `/usr/local/Q6/bin/`
- `/opt/Q6/bin/`
- `/cluster/software/Q6/bin/`
- `/home/username/software/Q6/bin/`

---

## Step 5: Add to settings.py

### Template for SLURM cluster

```python
CLUSTERS = {
    'localhost': {
        'Q_DIR': Q_DIR,
        'QDYN': os.path.join(Q_DIR, 'bin', 'qdyn'),
        'QPREP': os.path.join(Q_DIR, 'bin', 'qprep'),
        'QFEP': os.path.join(Q_DIR, 'bin', 'qfep'),
        'QCALC': os.path.join(Q_DIR, 'bin', 'qcalc'),
        'SUBMIT_CMD': 'bash',
        'MODULES': '',
        'NODES': 1,
        'NTASKS': 1,
        'CPUS_PER_TASK': 1,
        'TIME': '10:00:00',
        'PARTITION': 'normal',
    },
    'your_cluster': {
        'Q_DIR': '/cluster/software/Q6',
        'QDYN': '/cluster/software/Q6/bin/qdyn',
        'QPREP': '/cluster/software/Q6/bin/qprep',
        'QFEP': '/cluster/software/Q6/bin/qfep',
        'QCALC': '/cluster/software/Q6/bin/qcalc',
        'SUBMIT_CMD': 'sbatch',               # SLURM
        'MODULES': 'module load intel/2021 openmpi/4.0',
        'NODES': 1,
        'NTASKS': 8,                         # 8 parallel tasks
        'CPUS_PER_TASK': 2,                  # 2 CPUs per task → 16 total
        'TIME': '24:00:00',                  # HH:MM:SS format
        'PARTITION': 'gpu',                  # Specific partition if needed
    },
}
```

### Common Scheduler Fields

| Field | SLURM | PBS | SGE | Note |
|-------|-------|-----|-----|------|
| Submit command | `sbatch` | `qsub` | `qsub` | Used to launch job scripts |
| Parallel tasks | `#SBATCH -n` | `#PBS -l nodes=` | `#$ -pe` | Distributed computing |
| CPU count | `#SBATCH -c` | N/A | N/A | CPUs per task |
| Memory | `#SBATCH --mem` | `#PBS -l mem=` | `#$ -l h_vmem=` | Per node or per slot |
| Time limit | `#SBATCH -t` | `#PBS -l walltime=` | `#$ -l h_rt=` | HH:MM:SS format |
| Partition | `#SBATCH -p` | `#PBS -q` | N/A | Resource queue |

---

## Step 6: Create & Test Job Template

The submission script templates in `INPUTS/FEP_submit*.sh` will be auto-populated.

Test locally first:

```bash
# Example: Run simple Q command on cluster via settings.py
python -c "
import settings as s
cluster = s.CLUSTERS['your_cluster']
print('Q_DIR:', cluster['Q_DIR'])
print('QDYN:', cluster['QDYN'])
print('SUBMIT_CMD:', cluster['SUBMIT_CMD'])
"
```

Test job submission (dry-run):

```bash
# Generate a test FEP setup
python QligFEP.py -l1 lig1.pdb -l2 lig2.pdb -f OPLS2015 -S prot.pdb -C your_cluster -n 1 --no-submit

# Inspect generated job script
cat 17_22_FEP/FEP_submit.sh

# Verify script syntax
bash -n 17_22_FEP/FEP_submit.sh  # Should report no errors
```

Submit test job:

```bash
# Small test run
cd 17_22_FEP
sbatch FEP_submit.sh  # or qsub, as needed

# Monitor job
squeue -u username   # SLURM
qstat                 # PBS
```

---

## Step 7: Verify Job Execution

```bash
# Check job status
squeue -u username -j JOB_ID

# Check job output
cat slurm-JOB_ID.out
tail -f 17_22_FEP/md_1/.output_$name

# If job failed, check error logs
cat slurm-JOB_ID.err
ls -la 17_22_FEP/md_*/  # Check for output files
```

---

## Common Issues & Solutions

| Problem | Cause | Solution |
|---------|-------|----------|
| "sbatch: command not found" | SLURM not in PATH or not loaded | Check modules: `module load slurm` (if available) or verify cluster type |
| "Q: command not found" in job output | Modules not loaded in job script | Add module load commands to MODULES field or verify Q_DIR is correct |
| Job timeout | TIME limit too short | Increase TIME in cluster config (HH:MM:SS format) for larger runs |
| "Resource limit exceeded" | Asking for too many CPUs/memory | Reduce NTASKS or disable hyperthreading (CPUS_PER_TASK = 1) |
| Job stuck in queue | Partition unavailable or queue full | Check `sinfo` (SLURM) or `qstat` (PBS); try different PARTITION |
| Wrong Q executable version | Multiple Q installations on cluster | Verify full path in 'QDYN'/'QPREP'/'QFEP'/'QCALC' fields |

---

## After Configuration

1. **Update DEFAULT** if this is your main cluster:
   ```python
   DEFAULT = 'your_cluster'
   ```

2. **Test with small runs first**:
   ```bash
   python QligFEP.py ... -C your_cluster -n 1  # 1 lambda window
   python QresFEP.py ... -C your_cluster -n 1
   ```

3. **Monitor first job to completion**:
   ```bash
   squeue -u username
   tail -f output.log
   ```

4. **Document cluster quirks** (e.g., "Need to run from /cluster/scratch, not $HOME"):
   - Add as comment in settings.py
   - Share with team in README or wiki

---

## Reference

- **SLURM Documentation**: https://slurm.schedmd.com/sbatch.html
- **PBS Documentation**: https://www.altair.com/pbs-professional
- **SGE Documentation**: https://www.univa.com/resources/files/univa_grid_engine_overview_2019.pdf
- **settings.py Format**: See `.github/instructions/settings-py.instructions.md`

