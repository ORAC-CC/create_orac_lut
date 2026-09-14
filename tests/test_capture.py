from pathlib import Path

import pytest

from oraclut.validation.capture import (
    ReferenceNotUpdated,
    capture_verified_reference,
    fingerprint,
)


def test_capture_rejects_unchanged_source(tmp_path: Path):
    source = tmp_path / "source.nc"
    source.write_bytes(b"legacy")
    before = fingerprint(source)
    with pytest.raises(ReferenceNotUpdated):
        capture_verified_reference(
            source,
            tmp_path / "capture.nc",
            before,
            tmp_path / "manifest.json",
            configuration={"case": "test"},
            idl_version="IDL 8.9.0",
            mie_dlm_path="mie",
            disort_dlm_path="disort2",
        )
    assert not (tmp_path / "capture.nc").exists()


def test_capture_copies_only_an_updated_source_and_records_hash(tmp_path: Path):
    source = tmp_path / "source.nc"
    source.write_bytes(b"before")
    before = fingerprint(source)
    source.write_bytes(b"after")
    destination = tmp_path / "nested" / "capture.nc"
    manifest = tmp_path / "nested" / "manifest.json"
    record = capture_verified_reference(
        source,
        destination,
        before,
        manifest,
        configuration={"dimensions": {"channels": 1}},
        idl_version="IDL 8.9.0",
        mie_dlm_path="mie",
        disort_dlm_path="disort2",
    )
    assert destination.read_bytes() == b"after"
    assert record["source_size_bytes"] == len(b"after")
    assert record["source_sha256"] == fingerprint(destination).sha256
    assert manifest.exists()
