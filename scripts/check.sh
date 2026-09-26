#!/bin/sh
# Run from any directory; stop at the first failed gate.
set -eu
cd "$(dirname "$0")/.."

PYTHON=${PYTHON:-python3}
for verifier in K8 Sigma directional invisibility; do
    "$PYTHON" "paper/verify_${verifier}.py"
done
"$PYTHON" tests/regression_checks.py
"$PYTHON" threadB/reproduce_theorem4.py
sha256sum -c hashes.txt
