#!/usr/bin/env bash
set -euo pipefail

# Local repo hygiene helper:
# - reports large tracked files (current tree)
# - reports large files in git history
# - optional cleanup of local build artifacts

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THRESHOLD_MB="${THRESHOLD_MB:-50}"
THRESHOLD_BYTES=$((THRESHOLD_MB * 1024 * 1024))
DO_CLEAN=0

if [[ "${1:-}" == "--clean-local" ]]; then
  DO_CLEAN=1
fi

cd "${ROOT_DIR}"

echo "Repo: ${ROOT_DIR}"
echo "Threshold: ${THRESHOLD_MB} MB"
echo

echo "Largest tracked files in current tree:"
git ls-files | while read -r f; do
  [[ -f "${f}" ]] || continue
  size="$(stat -f '%z' "${f}")"
  if [[ "${size}" -ge "${THRESHOLD_BYTES}" ]]; then
    echo "${size}  ${f}"
  fi
done | sort -nr || true
echo

echo "Top files in git history:"
git rev-list --objects --all \
  | git cat-file --batch-check='%(objecttype) %(objectname) %(objectsize) %(rest)' \
  | awk '$1=="blob"{print $3 "\t" $4}' \
  | sort -nr \
  | head -n 20 || true
echo

if [[ "${DO_CLEAN}" -eq 1 ]]; then
  echo "Cleaning local artifacts (dist, iconset, pycache)..."
  rm -rf dist assets/InSituCore.iconset __pycache__ app/__pycache__ utils/__pycache__ utils/karospace/__pycache__ utils/mana/__pycache__
  echo "Done."
fi
