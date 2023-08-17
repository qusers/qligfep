import os
import subprocess
import sys

def run_analyze_FEP(parent_folder, command_args):
    for root, dirs, files in os.walk(parent_folder):
        for dir_name in dirs:
            if dir_name.startswith('FEP_'):
                full_path = os.path.join(root, dir_name)
                cmd = ['python', 'analyze_FEP.py', '-F', full_path] + command_args
                subprocess.run(cmd)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script_name.py <parent_folder> <command_args>")
        sys.exit(1)

    parent_folder = sys.argv[1]
    command_args = sys.argv[2:]

    if not os.path.isdir(parent_folder):
        print(f"Error: '{parent_folder}' is not a valid directory.")
        sys.exit(1)

    run_analyze_FEP(parent_folder, command_args)
