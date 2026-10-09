#!/usr/bin/env bash
# Create the annotated Git tag of a LUT code release (LUT_CODE_VERSIONS.md, section 5).
#
#   scripts/tag_lut_code_release.sh [-m MESSAGE_FILE] [--push]
#
# The tag name is taken from CODE_VERSION (code_release).  The script refuses
# to tag unless the working tree is clean (git status, including untracked
# files), HEAD is on main, CODE_VERSION parses (python -m oraclut.version),
# and the tag does not already exist locally or on origin.  With --push the
# commit and the tag are pushed to origin (never with --force).  The message
# file should describe the scientific content of the release and the product
# version it produces; without -m the editor opens.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$REPO_ROOT"
PYTHON="${ORACLUT_PYTHON:-/home/g/grainger/miniforge3/envs/science/bin/python}"
GIT="${ORACLUT_GIT:-git}"

message_file=""
push=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        -m) message_file="${2:?message file}"; shift 2 ;;
        --push) push=1; shift ;;
        -h|--help) sed -n '2,13p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "tag_lut_code_release.sh: unknown argument $1" >&2; exit 2 ;;
    esac
done

die() { echo "tag_lut_code_release.sh: $*" >&2; exit 1; }

tag="$(sed -n 's/^code_release *= *//p' CODE_VERSION | head -1)"
[[ -n "$tag" ]] || die "CODE_VERSION does not define code_release"
[[ "$tag" =~ ^lut-code-v[0-9]+(\.[0-9]+)*$ ]] || die "code_release '$tag' is not a lut-code-v<n>[.<m>] tag name"
PYTHONPATH=src "$PYTHON" -m oraclut.version > /dev/null || die "CODE_VERSION is not valid"

branch="$("$GIT" rev-parse --abbrev-ref HEAD)"
[[ "$branch" == "main" ]] || die "HEAD is on '$branch', not main"
# Clean means what the version banner means: no tracked change anywhere and no
# untracked file under the production paths (oraclut.version.PRODUCTION_PATHS).
status="$("$GIT" status --porcelain --untracked-files=no)"
[[ -z "$status" ]] || die "working tree has uncommitted tracked changes:"$'\n'"$status"
tree_state="$(PYTHONPATH=src "$PYTHON" -m oraclut.version --json | "$PYTHON" -c 'import json,sys; print(json.load(sys.stdin)["git"]["working_tree"])')"
[[ "$tree_state" == "CLEAN" ]] || die "working tree is $tree_state (see python -m oraclut.version)"
"$GIT" rev-parse --verify --quiet "refs/tags/$tag" > /dev/null && die "tag $tag already exists locally"
if "$GIT" ls-remote --exit-code --tags origin "refs/tags/$tag" > /dev/null 2>&1; then
    die "tag $tag already exists on origin"
fi

if [[ -n "$message_file" ]]; then
    [[ -s "$message_file" ]] || die "message file $message_file is empty"
    "$GIT" tag -a "$tag" -F "$message_file"
else
    "$GIT" tag -a "$tag"
fi
echo "created annotated tag $tag at $("$GIT" rev-parse HEAD)"
PYTHONPATH=src "$PYTHON" -m oraclut.version | sed -n '1,8p'

if [[ "$push" -eq 1 ]]; then
    "$GIT" push origin main
    "$GIT" push origin "$tag"
    "$GIT" ls-remote origin refs/heads/main "refs/tags/$tag"
fi
