# Installation Robustness Improvements

## Summary of Changes

We've made the QligFEP installation process more robust and portable across different machines. Here are the key improvements:

### 1. **Template-Based Configuration System**

**Before**: settings.py was edited directly with hardcoded paths
**After**: Uses a template file (`settings.py.template`) with placeholders that get substituted

**Files Created/Modified**:
- `settings.py.template` - Template with `$Q_PATH`, `$SCHROD_DIR`, `$CLUSTER_NAME` placeholders
- `settings.py` - Generated from template on each setup run
- Both scripts now restore from template before substituting values

### 2. **Python Setup Script** (`setup_qligfep.py`)

**New interactive setup tool with:**
- Auto-detection of local Q installation in repository
- Path validation (checks Q binaries actually exist)
- Error handling and user guidance
- Non-interactive mode for automation/CI
- Backwards compatible placeholder handling

**Usage**:
```bash
python setup_qligfep.py                                    # Interactive
python setup_qligfep.py /path/to/Q "" default_cluster   # Non-interactive
```

### 3. **Improved Bash Setup Script** (`qligfep_init.sh`)

**Enhancements**:
- Now uses template file instead of git checkout
- Updated placeholder handling
- Clearer prompts and validation
- Non-interactive mode support
- Works for environments without git

**Usage**:
```bash
bash qligfep_init.sh                              # Interactive
bash qligfep_init.sh /path/to/Q "" cluster_name # Non-interactive
```

### 4. **Installation Documentation** (`INSTALLATION.md`)

Comprehensive guide covering:
- Quick start with both setup scripts
- How the system works
- Configuration details (Q path, Schrodinger, clusters)
- Troubleshooting
- Multi-machine support

## Key Benefits

✅ **Portable**: Works on different machines without modification
✅ **Automatic**: Detects bundled Q installation automatically  
✅ **Flexible**: Supports custom Q paths, Schrodinger optional, multiple clusters
✅ **Safe**: Can run setup multiple times, template restores each time
✅ **Robust**: Both Python and Bash options, fallback paths
✅ **Validated**: Checks Q binaries exist before using
✅ **Documented**: Clear setup guide with examples

## For Users on Different Machines

### Local Q Installation (in repository):
```bash
python setup_qligfep.py
# Accepts auto-detected local Q path
```

### Q Installed Elsewhere:
```bash
python setup_qligfep.py /usr/local/q/bin "" default
```

### Multiple Clusters:
Edit `settings.py` after setup to add additional cluster configurations.

## Files Updated

| File | Purpose |
|------|---------|
| `settings.py.template` | Template with placeholders (NEW) |
| `setup_qligfep.py` | Python setup script (NEW) |
| `qligfep_init.sh` | Updated bash script |
| `INSTALLATION.md` | Setup documentation (NEW) |

## Technical Details

### Placeholder Substitution

The template uses clear placeholder names:
- `$Q_PATH` - Path to Q binaries directory
- `$SCHROD_DIR` - Schrödinger installation path
- `$CLUSTER_NAME` - Default HPC cluster name

Both setup scripts handle legacy placeholders for backwards compatibility.

### Template Restoration

On each setup run:
1. Copy `settings.py.template` to `settings.py`
2. Replace placeholders with actual values
3. Safe to run multiple times

### Validation

Python script validates:
- Q directory exists
- Q binaries (qprep/qdyn) are present
- Provides clear error messages if validation fails

## Testing

Current environment verified:
- ✅ Python setup script works
- ✅ Settings properly configured
- ✅ QresFEP runs without import errors
- ✅ Template system functional
