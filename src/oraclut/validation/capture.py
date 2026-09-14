"""Safe capture of a verified legacy NetCDF reference file."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil


@dataclass(frozen=True)
class FileFingerprint:
    path: str
    exists: bool
    size_bytes: int | None
    mtime_ns: int | None
    sha256: str | None


class ReferenceNotUpdated(RuntimeError):
    """Raised when the legacy source file was not changed by the run."""


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(path: str | Path) -> FileFingerprint:
    path = Path(path)
    if not path.is_file():
        return FileFingerprint(str(path), False, None, None, None)
    stat = path.stat()
    return FileFingerprint(str(path), True, stat.st_size, stat.st_mtime_ns, sha256_file(path))


def changed_since(before: FileFingerprint, after: FileFingerprint) -> bool:
    if not after.exists:
        return False
    if not before.exists:
        return True
    return (
        before.size_bytes != after.size_bytes
        or before.mtime_ns != after.mtime_ns
        or before.sha256 != after.sha256
    )


def capture_verified_reference(
    source: str | Path,
    destination: str | Path,
    before: FileFingerprint,
    manifest_path: str | Path,
    *,
    configuration: dict[str, object],
    idl_version: str,
    mie_dlm_path: str,
    disort_dlm_path: str,
) -> dict[str, object]:
    source = Path(source)
    destination = Path(destination)
    after = fingerprint(source)
    if not changed_since(before, after):
        raise ReferenceNotUpdated(
            f"Legacy source was not updated during this run: {source}"
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".tmp")
    shutil.copyfile(source, temporary)
    copied_hash = sha256_file(temporary)
    if copied_hash != after.sha256:
        temporary.unlink(missing_ok=True)
        raise IOError("Captured NetCDF hash does not match the legacy source")
    temporary.replace(destination)

    record: dict[str, object] = {
        "source_path": str(source),
        "captured_path": str(destination),
        "source_size_bytes": after.size_bytes,
        "source_sha256": after.sha256,
        "source_mtime_ns": after.mtime_ns,
        "capture_time_utc": datetime.now(timezone.utc).isoformat(),
        "idl_version": idl_version,
        "production_mie_dlm": mie_dlm_path,
        "production_disort_dlm": disort_dlm_path,
        "configuration": configuration,
        "before_fingerprint": asdict(before),
        "after_fingerprint": asdict(after),
    }
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(record, indent=2) + "\n")
    return record


def _load_fingerprint(path: str | Path) -> FileFingerprint:
    return FileFingerprint(**json.loads(Path(path).read_text()))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    fingerprint_parser = subparsers.add_parser("fingerprint")
    fingerprint_parser.add_argument("path", type=Path)
    fingerprint_parser.add_argument("--output", type=Path, required=True)

    capture_parser = subparsers.add_parser("capture")
    capture_parser.add_argument("source", type=Path)
    capture_parser.add_argument("destination", type=Path)
    capture_parser.add_argument("--before", type=Path, required=True)
    capture_parser.add_argument("--manifest", type=Path, required=True)
    capture_parser.add_argument("--configuration", type=Path, required=True)
    capture_parser.add_argument("--idl-version", required=True)
    capture_parser.add_argument("--mie-dlm-path", required=True)
    capture_parser.add_argument("--disort-dlm-path", required=True)
    args = parser.parse_args()

    if args.command == "fingerprint":
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(asdict(fingerprint(args.path)), indent=2) + "\n")
        return 0

    capture_verified_reference(
        args.source,
        args.destination,
        _load_fingerprint(args.before),
        args.manifest,
        configuration=json.loads(args.configuration.read_text()),
        idl_version=args.idl_version,
        mie_dlm_path=args.mie_dlm_path,
        disort_dlm_path=args.disort_dlm_path,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
