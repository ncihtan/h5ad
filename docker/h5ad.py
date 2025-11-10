# This script is used in the Docker image to run the validator.
# The main purpose is so that when validation or internal errors occur an exit code of 1 is returned.
# The Docker image is used by DCQC: https://github.com/Sage-Bionetworks-Workflows/py-dcqc

import subprocess
import sys

def run_h5ad():

    h5ad_path = sys.argv[1]

    result = subprocess.run(
        [
            'HTAN-h5ad-validate',
            h5ad_path,
            'output_file.txt'
        ],
        capture_output=True,
        text=True,
        check=True
    )

    print(result.stdout)

    failure_strings = ["Cellxgene run has errors.", "HTAN Validation Failed.", "An error occurred"]

    for string in failure_strings:
        if string in result.stdout:
            lines = []
            with open('output_file.txt', 'r', encoding="utf-8") as f:
                for line in f:
                    lines.append(line)
            print(lines)
            return 1

    return 0

if __name__ == '__main__':
    sys.exit(run_h5ad())
