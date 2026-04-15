# Installation & Setup Guide

## Quick Start

There are two ways to configure QligFEP on your machine:

### 1. Using Python Setup Script (Recommended)

The Python script provides better validation and interactive setup:

```bash
python setup_qligfep.py
```

This will:
- Auto-detect local Q installation in the repository
- Validate Q binaries are present
- Prompt for Schrodinger path (optional)
- Set default HPC cluster
- Display configuration summary

Or run in non-interactive mode with arguments:

```bash
python setup_qligfep.py /path/to/Q/bin /path/to/schrodinger default_cluster
```

### 2. Using Bash Setup Script

For shell-only environments:

```bash
bash qligfep_init.sh
```

Same features as the Python script, but in bash. If no arguments are provided, it will prompt interactively.

With arguments:
```bash
bash qligfep_init.sh /path/to/Q/bin /path/to/schrodinger cluster_name
```

## How It Works

Both scripts:

1. **Detect Local Q**: Check if Q is compiled in the repository (`/Q/bin/`), offer to use it
2. **Validate Q**: Ensure Q binaries (qprep, qdyn, etc.) are present
3. **Update Configuration**: Replace placeholders in `settings.py` with your paths
4. **Generate settings.py**: From the template file `settings.py.template`

### File Structure

- `settings.py.template` - Template with placeholders (source of truth)
- `settings.py` - Generated configuration (gitignored after first run)
- `setup_qligfep.py` - Python setup script
- `qligfep_init.sh` - Bash setup script

## Configuration Details

### Q Installation Paths

Q_PATH should point to the directory containing Q binaries. Examples:

- **Local repository**: `/path/to/qligfep/Q/bin`
- **System-wide**: `/usr/local/q/bin`
- **Custom location**: `/home/user/software/q6/bin`

### Schrodinger (Optional)

Required only for:
- OPLS parameter generation via `ffld_server`
- Protein preparation via `Protein Prep Wizard`

Leave empty if not installed or not needed.

### HPC Clusters

The default cluster name becomes the default in `settings.py`. To add more clusters, edit `settings.py` directly and add entries to the `Q_DIR` dictionary and cluster configuration dictionary.

Example adding a new cluster named "CLUSTER2":

```python
Q_DIR = {'default': '/path/to/q/bin/',
         'CLUSTER2': '/path/to/q/bin/',  # Can point to same Q or different
        }

CLUSTER2 = {'NODES': '2',
            'NTASKS': '16',
            ...
           }
```

Then use with: `python QligFEP.py ... -C CLUSTER2`

## Troubleshooting

### "Q binaries not found"

Ensure Q is properly installed/compiled at the specified path. Check for:
- `qprep` executable
- `qdyn` executable  
- `qfep` executable
- `qcalc` script/executable

###  Import errors after setup

If you still get import errors:

1. Make sure settings.py exists and is readable
2. Verify all paths in settings.py are correct
3. Try running setup again: `python setup_qligfep.py`

### Running multiple times

It's safe to run setup multiple times:
1. Template is restored each time
2. Placeholders are replaced with new values
3. Works with different Q paths for different environments

## For Different Machines

To install on a new machine:

1. Clone the repository
2. If Q is in the repo (in `/Q/bin/`), run setup and accept default
3. If Q is elsewhere, provide the path when prompted
4. Done! Ready to use

The template-based approach ensures users only need to run setup once, and `settings.py` is properly configured without breaking other installations.
