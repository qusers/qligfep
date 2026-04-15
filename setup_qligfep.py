#!/usr/bin/env python3
"""
QligFEP Configuration Setup Script

This script provides a robust way to configure QligFEP for your system,
detecting local Q installations and validating paths.
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Tuple


def find_q_binaries(search_dir: Path) -> bool:
    """Check if Q binaries (qprep, qdyn) exist in a directory."""
    return (search_dir / "qprep").exists() or (search_dir / "qdyn").exists()


def validate_q_path(q_path: str) -> Tuple[bool, str]:
    """Validate Q installation path."""
    q_dir = Path(q_path)
    
    if not q_dir.exists():
        return False, f"Directory does not exist: {q_path}"
    
    if not find_q_binaries(q_dir):
        return False, f"Q binaries (qprep/qdyn) not found in: {q_path}"
    
    return True, f"Q installation valid at: {q_path}"


def find_local_q_installation() -> Optional[Path]:
    """Attempt to find Q installed in the repository."""
    script_dir = Path(__file__).parent
    local_q_bin = script_dir / "Q" / "bin"
    
    if find_q_binaries(local_q_bin):
        return local_q_bin
    
    # Also check src/q6 for newly compiled binaries
    local_q_src = script_dir / "Q" / "src" / "q6"
    if find_q_binaries(local_q_src):
        return local_q_src
    
    return None


def restore_template_file(file_path: Path) -> bool:
    """Restore file from template if available."""
    template_path = file_path.parent / (file_path.name + ".template")
    
    if template_path.exists():
        try:
            shutil.copy(template_path, file_path)
            return True
        except Exception:
            return False
    
    return False


def update_settings_file(
    settings_path: Path,
    q_path: str,
    schrod_path: str,
    cluster_name: str
) -> bool:
    """Update settings.py with configuration values."""
    try:
        with open(settings_path, 'r') as f:
            content = f.read()
        
        # Replace placeholders (handle both old and new template formats)
        content = content.replace('$Q_PATH', q_path)
        content = content.replace('$PATH', q_path)  # For backwards compatibility
        content = content.replace('$SCHROD_DIR', schrod_path)
        content = content.replace('$SCHROD', schrod_path)  # For backwards compatibility
        content = content.replace('$CLUSTER_NAME', cluster_name)
        content = content.replace('$DEFAULT', cluster_name)  # For backwards compatibility
        
        with open(settings_path, 'w') as f:
            f.write(content)
        
        return True
    except Exception as e:
        print(f"Error updating settings.py: {e}")
        return False


def interactive_setup() -> Tuple[str, str, str]:
    """Run interactive setup prompts."""
    print("\n" + "="*50)
    print("QligFEP Configuration Setup")
    print("="*50 + "\n")
    
    # Q Path
    print("Step 1: Locate Q Installation")
    print("-" * 50)
    
    local_q = find_local_q_installation()
    if local_q:
        print(f"Found local Q installation: {local_q}")
        use_local = input("Use this Q installation? (y/n) [y]: ").strip().lower()
        if use_local != 'n':
            q_path = str(local_q)
        else:
            q_path = input("Enter absolute path to Q directory: ").strip()
    else:
        print("No local Q installation found.")
        q_path = input("Enter absolute path to Q directory: ").strip()
    
    # Validate Q path
    while True:
        valid, msg = validate_q_path(q_path)
        print(msg)
        if valid:
            break
        q_path = input("Enter absolute path to Q directory (or press Ctrl+C to exit): ").strip()
    
    # Schrodinger Path (optional)
    print("\nStep 2: Schrodinger Installation (optional)")
    print("-" * 50)
    schrod_path = input("Enter absolute path to Schrodinger (leave blank if not installed): ").strip()
    
    if schrod_path and not Path(schrod_path).exists():
        print(f"Warning: {schrod_path} does not exist. Proceeding anyway...")
    
    # Cluster Name
    print("\nStep 3: Default HPC Cluster")
    print("-" * 50)
    cluster_name = input("Enter default HPC cluster name [default]: ").strip() or "default"
    
    return q_path, schrod_path, cluster_name


def main():
    """Main setup function."""
    script_dir = Path(__file__).parent
    settings_file = script_dir / "settings.py"
    
    print("\n" + "="*50)
    print("QligFEP Installation Setup")
    print("="*50)
    
    if len(sys.argv) > 1:
        # Command-line mode
        if len(sys.argv) < 4:
            print("Usage: python setup_qligfep.py <Q_PATH> <SCHROD_PATH> <CLUSTER_NAME>")
            print("       or: python setup_qligfep.py  (for interactive mode)")
            sys.exit(1)
        
        q_path = sys.argv[1]
        schrod_path = sys.argv[2] if sys.argv[2] != "none" else ""
        cluster_name = sys.argv[3]
    else:
        # Interactive mode
        try:
            q_path, schrod_path, cluster_name = interactive_setup()
        except KeyboardInterrupt:
            print("\n\nSetup cancelled.")
            sys.exit(0)
    
    # Restore template if in git repo
    print("\nRestoring settings template...", end="")
    if restore_template_file(settings_file):
        print(" ✓")
    else:
        print(" (not a git repo, proceeding)")
    
    # Update settings
    print("Updating configuration...", end="")
    if update_settings_file(settings_file, q_path, schrod_path, cluster_name):
        print(" ✓")
    else:
        print(" ✗")
        print("Failed to update settings.py")
        sys.exit(1)
    
    # Display summary
    print("\n" + "="*50)
    print("Setup Complete!")
    print("="*50)
    print(f"Q directory:           {q_path}")
    print(f"Schrodinger directory: {schrod_path or '(not installed)'}")
    print(f"Default cluster:       {cluster_name}")
    print("\nYou can now run:")
    print("  python QligFEP.py -h")
    print("  python QresFEP.py -h")
    print("  python QLIE.py -h")
    print("\nTo modify HPC clusters, edit settings.py directly.")
    print()


if __name__ == "__main__":
    main()
