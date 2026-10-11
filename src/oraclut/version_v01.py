"""Code-version identification of the ORAC LUT generator.

One authoritative definition, the ``CODE_VERSION`` file at the repository
root, declares the source release (an annotated Git tag ``lut-code-v<n>``),
the LUT product version that release produces, any other product versions it
is explicitly allowed to produce, and a short identifier of the scientific
configuration.  Git supplies the exact commit, the relation of HEAD to the
declared release tag and whether the working tree is modified.  This module
reads both, prints the start-up banner, enforces the agreement between the
requested LUT version and the source release, and provides the command line

    python -m oraclut.version                 banner
    python -m oraclut.version --json          machine-readable record
    python -m oraclut.version --list-releases historical lut-code-v* releases
    python -m oraclut.version --verify PROVENANCE.json [--lut FILE]
                                              which revision generated an output

Nothing here touches a scientific calculation.  See LUT_CODE_VERSIONS.md.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import shutil
import socket
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

REPO_ROOT = Path(__file__).resolve().parents[2]
CODE_VERSION_FILE = REPO_ROOT / "CODE_VERSION"
RELEASE_TAG_PATTERN = re.compile(r"^lut-code-v\d+(\.\d+)*$")
# Untracked files under these paths make the working tree MODIFIED: they are
# the paths whose content defines a LUT product (scripts/submit_v25_cloud_luts.sh
# uses the same notion of production content).
PRODUCTION_PATHS = ("CODE_VERSION", "create_orac_luts.py", "makerunfile", "src", "mie",
                    "create_orac_lut", "runs", "configs", "scripts")
_GIT_TIMEOUT = 60


class VersionError(RuntimeError):
    """The version definition is missing, malformed or violated."""


class LutVersionMismatch(VersionError):
    """A LUT version was requested that this source release does not produce."""


# ---------------------------------------------------------------------------
# CODE_VERSION
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CodeVersion:
    """The declared source release (contents of ``CODE_VERSION``)."""

    code_release: str
    lut_version: int
    compatible_lut_versions: tuple[int, ...]
    scientific_config: str
    scientific_config_description: str
    source: str

    def accepts(self, lut_version: int | None) -> bool:
        """Whether this release may produce ``lut_version``."""

        return lut_version is not None and int(lut_version) in (self.lut_version, *self.compatible_lut_versions)


def read_code_version(path: str | os.PathLike = CODE_VERSION_FILE) -> CodeVersion:
    """Parse ``CODE_VERSION``; every field is required and validated."""

    path = Path(path)
    if not path.is_file():
        raise VersionError(f"version definition not found: {path}")
    values: dict[str, str] = {}
    for number, raw in enumerate(path.read_text().splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise VersionError(f"{path}:{number}: expected 'name = value', got {raw!r}")
        name, value = (part.strip() for part in line.split("=", 1))
        if name in values:
            raise VersionError(f"{path}:{number}: {name!r} is set twice")
        values[name] = value
    required = ("code_release", "lut_version", "compatible_lut_versions", "scientific_config",
                "scientific_config_description")
    missing = [name for name in required if name not in values]
    if missing:
        raise VersionError(f"{path}: missing settings: {', '.join(missing)}")
    unknown = sorted(set(values) - set(required))
    if unknown:
        raise VersionError(f"{path}: unknown settings: {', '.join(unknown)}")
    release = values["code_release"]
    if not RELEASE_TAG_PATTERN.match(release):
        raise VersionError(f"{path}: code_release must be an annotated tag name like lut-code-v25 or "
                           f"lut-code-v25.1, got {release!r}")
    try:
        lut_version = int(values["lut_version"])
        compatible = tuple(int(item) for item in re.split(r"[\s,]+", values["compatible_lut_versions"]) if item)
    except ValueError as exc:
        raise VersionError(f"{path}: lut_version and compatible_lut_versions must be integers") from exc
    if lut_version in compatible:
        raise VersionError(f"{path}: compatible_lut_versions must not repeat lut_version {lut_version}")
    if not values["scientific_config"]:
        raise VersionError(f"{path}: scientific_config must not be empty")
    return CodeVersion(release, lut_version, compatible, values["scientific_config"],
                       values["scientific_config_description"], str(path))


# ---------------------------------------------------------------------------
# Git metadata
# ---------------------------------------------------------------------------

def _git_executable() -> str:
    return os.environ.get("ORACLUT_GIT") or shutil.which("git") or "/usr/bin/git"


def _git(root: Path, *args: str) -> str:
    env = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
    env["LC_ALL"] = "C"
    completed = subprocess.run([_git_executable(), "-C", str(root), *args], capture_output=True, text=True,
                               timeout=_GIT_TIMEOUT, env=env)
    if completed.returncode != 0:
        raise VersionError(f"git {' '.join(args)}: {completed.stderr.strip() or completed.stdout.strip()}")
    return completed.stdout.rstrip("\n")


@dataclass(frozen=True)
class GitState:
    """What Git reports about the source tree a calculation runs from."""

    available: bool
    detail: str = ""                 # why Git metadata is unavailable
    commit: str | None = None        # full hash of HEAD
    branch: str | None = None        # branch name, or "detached"
    describe: str | None = None      # git describe --tags --long --always
    tracked_changes: tuple[str, ...] = ()
    untracked_production: tuple[str, ...] = ()
    release: str | None = None               # the release tag the two fields below refer to
    release_tag_commit: str | None = None    # commit that tag names, if the tag exists
    commits_after_release: int | None = None

    @property
    def short_commit(self) -> str:
        return self.commit[:7] if self.commit else "unavailable"

    @property
    def modified(self) -> bool | None:
        if not self.available:
            return None
        return bool(self.tracked_changes or self.untracked_production)

    @property
    def working_tree(self) -> str:
        if not self.available:
            return f"UNKNOWN (Git metadata unavailable: {self.detail})"
        return "MODIFIED" if self.modified else "CLEAN"

    def release_relation(self, release: str | None = None) -> str:
        """Relation of HEAD to the release tag this state was queried for."""

        release = release or self.release or "(no release declared)"
        if not self.available:
            return "unknown (Git metadata unavailable)"
        if self.release != release or self.release_tag_commit is None:
            return f"tag {release} not present in this repository"
        if self.release_tag_commit == self.commit:
            return f"HEAD is the tagged release {release}"
        if self.commits_after_release is not None:
            return f"HEAD is {self.commits_after_release} commit(s) after tag {release}"
        return f"HEAD is not the tagged release {release} ({self.release_tag_commit[:7]})"


def git_state(root: str | os.PathLike = REPO_ROOT, release: str | None = None) -> GitState:
    """Query Git for ``root``; never raises, reports unavailability instead."""

    root = Path(root)
    try:
        commit = _git(root, "rev-parse", "HEAD")
        branch = _git(root, "rev-parse", "--abbrev-ref", "HEAD")
        if branch == "HEAD":
            branch = "detached"
        try:
            describe = _git(root, "describe", "--tags", "--long", "--always", "--match", "lut-code-v*")
        except VersionError:
            describe = commit[:7]
        status = _git(root, "status", "--porcelain", "--untracked-files=all")
        tracked, untracked = [], []
        for line in status.splitlines():
            code, path = line[:2], line[3:]
            if " -> " in path:
                path = path.split(" -> ", 1)[1]
            if code == "??":
                if any(path == p or path.startswith(p + "/") for p in PRODUCTION_PATHS):
                    untracked.append(path)
            else:
                tracked.append(path)
        tag_commit = None
        after = None
        if release:
            try:
                tag_commit = _git(root, "rev-parse", "--verify", "--quiet", f"refs/tags/{release}^{{commit}}")
            except VersionError:
                tag_commit = None
            if tag_commit and tag_commit != commit:
                try:
                    after = int(_git(root, "rev-list", "--count", f"{tag_commit}..{commit}"))
                except VersionError:
                    after = None
        return GitState(True, "", commit, branch, describe, tuple(tracked), tuple(untracked), release, tag_commit, after)
    except (VersionError, OSError, subprocess.TimeoutExpired) as exc:
        return GitState(False, str(exc).splitlines()[0] if str(exc) else type(exc).__name__)


# ---------------------------------------------------------------------------
# Summary, banner, consistency check, provenance record
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class VersionSummary:
    code: CodeVersion
    git: GitState

    def banner(self, *, lut_version: int | None = None, generator: str | None = None, width: int = 60) -> str:
        rule = "-" * width
        lines = ["ORAC LUT Generator", rule]
        if generator:
            lines.append(f"Generator:         {generator}")
        requested = "" if lut_version is None else f" (requested {int(lut_version)})"
        lines += [
            f"LUT version:       {self.code.lut_version}{requested}",
            f"Code release:      {self.code.code_release}",
            f"Git commit:        {self.git.short_commit}"
            + (f" ({self.git.commit})" if self.git.commit else ""),
            f"Git describe:      {self.git.describe or 'unavailable'}",
            f"Release tag:       {self.git.release_relation(self.code.code_release)}",
            f"Working tree:      {self.git.working_tree}",
        ]
        for path in self.git.tracked_changes[:8]:
            lines.append(f"    modified:      {path}")
        for path in self.git.untracked_production[:8]:
            lines.append(f"    untracked:     {path}")
        hidden = len(self.git.tracked_changes) + len(self.git.untracked_production) - 16
        if hidden > 0:
            lines.append(f"    ... and {hidden} more")
        lines += [
            f"Scientific config: {self.code.scientific_config}",
            f"                   {self.code.scientific_config_description}",
            f"Compatible LUTs:   {' '.join(str(v) for v in self.code.compatible_lut_versions) or 'none'}",
            rule,
        ]
        return "\n".join(lines)

    def record(self) -> dict:
        """JSON-able record of the source identity and software environment."""

        return {
            "code_release": self.code.code_release,
            "lut_version": self.code.lut_version,
            "compatible_lut_versions": list(self.code.compatible_lut_versions),
            "scientific_config": self.code.scientific_config,
            "scientific_config_description": self.code.scientific_config_description,
            "git": {
                "available": self.git.available,
                "commit": self.git.commit,
                "branch": self.git.branch,
                "describe": self.git.describe,
                "working_tree": self.git.working_tree,
                "modified": self.git.modified,
                "tracked_changes": list(self.git.tracked_changes),
                "untracked_production_files": list(self.git.untracked_production),
                "release_tag_commit": self.git.release_tag_commit,
                "release_relation": self.git.release_relation(self.code.code_release),
                "detail": self.git.detail,
            },
            "software": software_versions(),
        }


def software_versions() -> dict:
    """Versions of the interpreter, libraries and Fortran compiler in use."""

    versions: dict = {"python": sys.version.split()[0], "platform": platform.platform(),
                      "hostname": socket.gethostname()}
    try:
        import numpy
        versions["numpy"] = numpy.__version__
    except Exception:   # pragma: no cover - environment dependent
        pass
    try:
        import netCDF4
        versions["netCDF4"] = netCDF4.__version__
        versions["libnetcdf"] = netCDF4.__netcdf4libversion__
        versions["hdf5"] = netCDF4.__hdf5libversion__
    except Exception:   # pragma: no cover - environment dependent
        pass
    compiler = os.environ.get("ORACLUT_FORTRAN", "gfortran")
    versions["fortran_compiler"] = compiler
    versions["fortran_flags"] = os.environ.get("ORACLUT_FORTRAN_FLAGS", "")
    try:
        completed = subprocess.run([compiler, "--version"], capture_output=True, text=True, timeout=20)
        versions["fortran_compiler_version"] = (completed.stdout or completed.stderr).splitlines()[0]
    except Exception:
        versions["fortran_compiler_version"] = "unavailable"
    return versions


def summarise(root: str | os.PathLike = REPO_ROOT) -> VersionSummary:
    root = Path(root)
    code = read_code_version(root / "CODE_VERSION")
    return VersionSummary(code, git_state(root, code.code_release))


def print_banner(*, lut_version: int | None = None, generator: str | None = None,
                 root: str | os.PathLike = REPO_ROOT, file=None) -> VersionSummary:
    """Print the start-up version summary and return it."""

    summary = summarise(root)
    print(summary.banner(lut_version=lut_version, generator=generator), file=file or sys.stdout, flush=True)
    return summary


def check_lut_version(requested: int | None, code: CodeVersion, *, source: str = "the configuration") -> int:
    """Return ``requested`` if this release may produce it; raise LutVersionMismatch otherwise."""

    if requested is None:
        raise LutVersionMismatch(f"{source} does not state the LUT version; this source release "
                                 f"{code.code_release} produces LUT version {code.lut_version}")
    requested = int(requested)
    if code.accepts(requested):
        return requested
    compatible = ", ".join(str(v) for v in code.compatible_lut_versions) or "none"
    raise LutVersionMismatch(
        f"{source} requests LUT version {requested}, but this source release {code.code_release} produces "
        f"LUT version {code.lut_version} (compatible: {compatible}).  A V{requested} product must be generated "
        f"from the release that defines V{requested} (python -m oraclut.version --list-releases; "
        f"git worktree add <dir> lut-code-v{requested}), or the requested version must be corrected.  "
        "The output filename is never adjusted to hide this.  See LUT_CODE_VERSIONS.md."
    )


# ---------------------------------------------------------------------------
# Releases and provenance verification
# ---------------------------------------------------------------------------

def list_releases(root: str | os.PathLike = REPO_ROOT) -> list[dict]:
    """Annotated ``lut-code-v*`` tags with their commits, dates and subjects."""

    out = _git(Path(root), "tag", "-l", "lut-code-v*", "--sort=version:refname",
               "--format=%(refname:strip=2)%09%(*objectname)%09%(objectname)%09%(taggerdate:short)%09%(contents:subject)")
    releases = []
    for line in out.splitlines():
        tag, commit, obj, date, subject = (line.split("\t") + [""] * 5)[:5]
        releases.append({"tag": tag, "commit": commit or obj, "annotated": bool(commit), "date": date,
                         "subject": subject})
    return releases


def sha256_file(path: str | os.PathLike, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_provenance(provenance_path: str | os.PathLike, *, lut_path: str | os.PathLike | None = None,
                      root: str | os.PathLike = REPO_ROOT) -> dict:
    """Relate a provenance record to this repository and (optionally) to the output file."""

    root = Path(root)
    record = json.loads(Path(provenance_path).read_text())
    source = record.get("source", {})
    commit = (source.get("git") or {}).get("commit")
    result: dict = {"provenance": str(provenance_path), "code_release": source.get("code_release"),
                    "lut_version": record.get("lut_version"), "commit": commit,
                    "working_tree": (source.get("git") or {}).get("working_tree"),
                    "generated_utc": record.get("generated_utc")}
    if commit:
        try:
            _git(root, "cat-file", "-e", f"{commit}^{{commit}}")
            result["commit_in_repository"] = True
            tags = _git(root, "tag", "--points-at", commit, "-l", "lut-code-v*").split()
            containing = _git(root, "tag", "--contains", commit, "-l", "lut-code-v*").split()
            result["release_tags_at_commit"] = tags
            result["release_tags_containing_commit"] = containing
            result["head_is_commit"] = _git(root, "rev-parse", "HEAD") == commit
        except VersionError as exc:
            result["commit_in_repository"] = False
            result["detail"] = str(exc)
    else:
        result["commit_in_repository"] = None
    output = record.get("output", {})
    if lut_path is None and output.get("file"):
        candidate = Path(provenance_path).with_name(output["file"])
        lut_path = candidate if candidate.exists() else None
    if lut_path is not None and Path(lut_path).exists():
        digest = sha256_file(lut_path)
        result["lut_file"] = str(lut_path)
        result["lut_sha256_matches"] = digest == output.get("sha256")
    result["consistent"] = bool(result.get("commit_in_repository")) and result.get("working_tree") == "CLEAN" \
        and result.get("lut_sha256_matches", True)
    return result


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m oraclut.version",
                                     description="Show the ORAC LUT generator's code version and Git state.")
    parser.add_argument("--root", type=Path, default=REPO_ROOT, help="repository root (default: this checkout)")
    parser.add_argument("--json", action="store_true", help="print the machine-readable record")
    parser.add_argument("--list-releases", action="store_true", help="list the lut-code-v* releases")
    parser.add_argument("--verify", type=Path, metavar="PROVENANCE.json",
                        help="relate a provenance record (and its LUT) to this repository")
    parser.add_argument("--lut", type=Path, help="LUT file to check against the provenance record")
    parser.add_argument("--check-lut-version", type=int, metavar="N",
                        help="exit 1 unless this release may produce LUT version N")
    args = parser.parse_args(argv)
    try:
        if args.list_releases:
            releases = list_releases(args.root)
            if not releases:
                print("no lut-code-v* releases found")
            for item in releases:
                print(f"{item['tag']:<18} {item['commit'][:12]}  {item['date']:<10}  {item['subject']}")
            return 0
        if args.verify:
            result = verify_provenance(args.verify, lut_path=args.lut, root=args.root)
            print(json.dumps(result, indent=1))
            return 0 if result["consistent"] else 1
        summary = summarise(args.root)
        if args.json:
            print(json.dumps(summary.record(), indent=1))
        else:
            print(summary.banner())
        if args.check_lut_version is not None:
            check_lut_version(args.check_lut_version, summary.code, source="the command line")
        return 0
    except VersionError as exc:
        print(f"oraclut.version: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
