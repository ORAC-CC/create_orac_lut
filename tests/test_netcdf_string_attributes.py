"""NetCDF text-attribute contract of the ORAC LUT writer and the metadata-only repair.

ORAC reads the axis ``spacing`` attributes with ``nc_get_att_string``, which
accepts NC_STRING only.  These tests check, with the NetCDF C library rather
than the netCDF4 package's decoded view, that the production writer emits
NC_STRING for every text attribute, that a plain ``setncattr`` reproduces the
production failure, and that the repair utility converts existing files while
preserving every byte of data and every unrelated piece of metadata.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import time
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

from oraclut.io import netcdf_c as nc
from oraclut.io import string_attributes as sa
from oraclut.io.v2 import put_attribute, write_v2_lut

V11_REFERENCE = Path("/network/group/aopp/eodg/RGG004_GRAINGER_ORACFILE/ORAC_LUTS/"
                     "aqua_modis_m_liquid-water_a01_p240_v11.nc")
OLD = time.time() - 24 * 3600


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def text_types(path: Path) -> dict[str, str]:
    return {a.location: a.type_name for a in nc.list_attributes(path)
            if a.xtype in (nc.NC_CHAR, nc.NC_STRING) and a.name != "_FillValue"}


def attribute_order(path: Path) -> list[tuple[str, str]]:
    return [(a.variable or "GLOBAL", a.name) for a in nc.list_attributes(path)]


def nan_with_payload() -> np.float32:
    bits = np.array([0x7FC12345], dtype=np.uint32)
    return bits.view(np.float32)[0]


def write_old_style_lut(path: Path) -> None:
    """A small LUT written the way the V22-V25 generator wrote it (text as NC_CHAR).

    It mixes everything the repair must cope with: interleaved numeric
    attributes, global text and numeric attributes, an attribute that is
    already NC_STRING, a character variable, a fill value, a NaN payload,
    a signed zero, and a chunked, compressed, shuffled variable.
    """

    with Dataset(path, "w", format="NETCDF4") as ds:
        for name, size in (("optical_depth", 3), ("effective_radius", 2), ("satellite_zenith", 2),
                           ("solar_zenith", 2), ("relative_azimuth", 2), ("channels", 2), ("st2", 4)):
            ds.createDimension(name, size)
        axes = {"optical_depth": ("uneven_logarithmic", [0.1, 1.0, 10.0]), "effective_radius": ("uneven_linear", [1, 2]),
                "satellite_zenith": ("linear", [0, 10]), "solar_zenith": ("linear", [0, 10]),
                "relative_azimuth": ("linear", [0, 180])}
        for name, (spacing, values) in axes.items():
            v = ds.createVariable(name, "f4", (name,))
            v[:] = values
            v.setncattr("long_name", name.replace("_", " "))
            v.setncattr("spacing", spacing)
            v.setncattr("units", "dimensionless")
            v.setncattr("valid_range", np.array([0.0, 1000.0], dtype=np.float32))
        platform = ds.createVariable("platform", "S1", ("st2",))
        platform[:] = np.array(list(b"aqua"), dtype="S1")
        platform.setncattr("long_name", "platform name")
        scalar = ds.createVariable("max_sat_zenith", "f4")
        scalar[...] = 90.0
        scalar.setncattr("units", "degrees")
        scalar.setncattr("valid_range", np.array([0.0, 90.0], dtype=np.float32))
        scalar.setncattr_string("already_string", "kept")
        scalar.setncattr("comment", "after the string one")
        t_dd = ds.createVariable("T_dd", "f4", ("optical_depth", "channels"), fill_value=np.float32(-999.0))
        t_dd[:] = np.array([[1.0, -0.0], [nan_with_payload(), 2.0], [-999.0, 3.0]], dtype=np.float32)
        t_dd.setncattr("long_name", "diffuse transmission")
        t_dd.setncattr("units", "dimensionless")
        r_0v = ds.createVariable("R_0v", "f4", ("relative_azimuth", "satellite_zenith", "solar_zenith",
                                                 "optical_depth", "effective_radius", "channels"),
                                 zlib=True, complevel=4, shuffle=True, chunksizes=(1, 2, 2, 3, 2, 2))
        r_0v[:] = np.arange(r_0v.size, dtype=np.float32).reshape(r_0v.shape) / 7.0
        r_0v.setncattr("long_name", "bi-directional reflectance")
        r_0v.setncattr("units", "dimensionless")
        r_0v.setncattr("valid_range", np.array([0.0, 1.0], dtype=np.float32))
        ds.setncattr("cloud_vertical_profile", "wet_adiabat")
        ds.setncattr("cloud_top_temperature_K", np.float32(240.0))
        ds.setncattr("cloud_emission_layers", np.int32(99))
        ds.setncattr("cloud_vertical_profile_note", "n" * 605)
    os.utime(path, (OLD, OLD))


@pytest.fixture
def old_lut(tmp_path):
    path = tmp_path / "aqua_modis_m_liquid-water_a01_p240_v25.nc"
    write_old_style_lut(path)
    reference = tmp_path / "reference_copy.nc"
    shutil.copyfile(path, reference)
    os.utime(reference, (OLD, OLD))
    return path, reference


def write_new_style_lut(path: Path) -> None:
    write_v2_lut(
        path, lut_level=2, revision=25,
        dimensions={"optical_depth": 2, "effective_radius": 2, "channels": 1},
        variables={"optical_depth": np.array([1.0, 2.0], "f4"), "effective_radius": np.array([5.0, 10.0], "f4"),
                   "T_dd": np.ones((2, 1), "f4")},
        variable_dimensions={"optical_depth": ("optical_depth",), "effective_radius": ("effective_radius",),
                             "T_dd": ("optical_depth", "channels")},
        variable_attributes={
            "optical_depth": {"long_name": "optical depth", "spacing": "logarithmic", "units": "dimensionless",
                              "valid_range": np.array([0.0, 100.0], "f4")},
            "effective_radius": {"long_name": "particle effective radius", "spacing": "linear", "units": "microns"},
            "T_dd": {"long_name": "diffuse transmission", "units": "dimensionless"}},
        global_attributes={"cloud_vertical_profile": "cirrostratus", "cloud_emission_layers": np.int32(99),
                           "cloud_vertical_profile_note": "x" * 600})


# ---------------------------------------------------------------------------
# The writer
# ---------------------------------------------------------------------------

def test_writer_emits_nc_string_for_every_text_attribute(tmp_path):
    path = tmp_path / "new.nc"
    write_new_style_lut(path)
    types = text_types(path)
    assert types and set(types.values()) == {"NC_STRING"}
    assert types["optical_depth:spacing"] == "NC_STRING"
    assert types["effective_radius:spacing"] == "NC_STRING"
    assert types["GLOBAL:cloud_vertical_profile"] == "NC_STRING"
    numeric = {a.location: a.type_name for a in nc.list_attributes(path) if a.xtype not in (nc.NC_CHAR, nc.NC_STRING)}
    assert numeric == {"optical_depth:valid_range": "NC_FLOAT", "GLOBAL:cloud_emission_layers": "NC_INT"}
    # The attribute ORAC reads, read the way ORAC reads it.
    assert nc.read_string_attribute_like_orac(path, "optical_depth", "spacing") == "logarithmic"
    assert nc.read_string_attribute_like_orac(path, "effective_radius", "spacing") == "linear"
    assert sa.audit_file(path).compliant


def test_plain_setncattr_reproduces_the_production_failure(tmp_path):
    path = tmp_path / "old.nc"
    write_old_style_lut(path)
    assert text_types(path)["optical_depth:spacing"] == "NC_CHAR"
    with pytest.raises(nc.NetCDFError) as info:
        nc.read_string_attribute_like_orac(path, "optical_depth", "spacing")
    assert info.value.status == nc.NC_ECHAR
    assert "Attempt to convert between text & numbers" in str(info.value)


def test_writer_fails_loudly_if_text_attributes_become_nc_char(tmp_path, monkeypatch):
    """The self-check must fail the write if a library default ever yields NC_CHAR again."""

    import oraclut.io.v2 as v2

    monkeypatch.setattr(v2, "put_attribute", lambda target, name, value: target.setncattr(name, value))
    with pytest.raises(ValueError, match="NC_STRING contract.*spacing is NC_CHAR"):
        write_new_style_lut(tmp_path / "regressed.nc")


def test_put_attribute_types_are_explicit(tmp_path):
    path = tmp_path / "attrs.nc"
    with Dataset(path, "w", format="NETCDF4") as ds:
        put_attribute(ds, "single", "text")
        put_attribute(ds, "several", ["a", "b"])
        put_attribute(ds, "numpy_text", np.array(["c", "d"]))
        put_attribute(ds, "number", np.float32(1.5))
        put_attribute(ds, "numbers", np.array([1, 2], dtype=np.int16))
        with pytest.raises(TypeError, match="bytes"):
            put_attribute(ds, "raw", b"bytes")
    types = {a.name: (a.type_name, a.length) for a in nc.list_attributes(path)}
    assert types == {"single": ("NC_STRING", 1), "several": ("NC_STRING", 2), "numpy_text": ("NC_STRING", 2),
                     "number": ("NC_FLOAT", 1), "numbers": ("NC_SHORT", 2)}


# ---------------------------------------------------------------------------
# The repair
# ---------------------------------------------------------------------------

def test_audit_lists_conversions_and_the_orac_requirement(old_lut):
    path, _ = old_lut
    audit = sa.audit_file(path)
    assert audit.file_format == "NC_FORMAT_NETCDF4"
    assert not audit.problems
    assert {c.location for c in audit.required_conversions} == {f"{axis}:spacing" for axis in sa.ORAC_AXIS_VARIABLES}
    assert "GLOBAL:cloud_vertical_profile_note" in {c.location for c in audit.conversions}
    assert "max_sat_zenith:already_string" not in {c.location for c in audit.conversions}
    assert "T_dd:_FillValue" not in {c.location for c in audit.conversions}
    assert audit.text_attribute_types == {"NC_CHAR": len(audit.conversions), "NC_STRING": 1}


def test_repair_converts_text_attributes_and_preserves_everything_else(old_lut):
    path, reference = old_lut
    order_before = attribute_order(path)
    with Dataset(path) as ds:
        storage_before = {n: (v.chunking(), v.filters(), v.endian(), v.dtype, v.shape) for n, v in ds.variables.items()}
    result = sa.repair_file(path, min_age_seconds=600)
    assert result.status == "repaired", result.message
    assert result.required_conversions == 5
    assert set(text_types(path).values()) == {"NC_STRING"}
    assert attribute_order(path) == order_before
    with Dataset(path) as ds:
        assert {n: (v.chunking(), v.filters(), v.endian(), v.dtype, v.shape) for n, v in ds.variables.items()} == storage_before
        t_dd = ds.variables["T_dd"]
        t_dd.set_auto_maskandscale(False)
        words = t_dd[:].view(np.uint32).ravel().tolist()
        assert words == [0x3F800000, 0x80000000, 0x7FC12345, 0x40000000, 0xC479C000, 0x40400000]
        assert ds.variables["T_dd"]._FillValue == np.float32(-999.0)
        assert ds.variables["max_sat_zenith"].getncattr("already_string") == "kept"
        assert ds.getncattr("cloud_emission_layers") == 99 and ds.getncattr("cloud_emission_layers").dtype == np.int32
    for axis in sa.ORAC_AXIS_VARIABLES:
        assert nc.read_string_attribute_like_orac(path, axis, "spacing") == result.orac_spacing_read[axis]
    assert result.orac_spacing_read["optical_depth"] == "uneven_logarithmic"
    audit = sa.audit_file(reference)
    report = sa.verify_preserved(reference, path, audit.conversions)
    assert report.ok, report.differences
    assert report.attributes_converted == len(audit.conversions)
    assert report.bytes_compared > 0
    assert not path.with_name(path.name + sa.TEMP_SUFFIX).exists()


def test_repair_is_idempotent(old_lut):
    path, _ = old_lut
    assert sa.repair_file(path, min_age_seconds=600).status == "repaired"
    digest, mtime = sha256(path), path.stat().st_mtime
    second = sa.repair_file(path, min_age_seconds=0)
    assert second.status == "compliant"
    assert sha256(path) == digest and path.stat().st_mtime == mtime
    assert second.orac_spacing_read["relative_azimuth"] == "linear"


def test_already_correct_file_is_reported_compliant_and_untouched(tmp_path):
    path = tmp_path / "new.nc"
    write_new_style_lut(path)
    digest = sha256(path)
    result = sa.repair_file(path, min_age_seconds=0)
    assert result.status == "compliant" and not result.conversions
    assert sha256(path) == digest


def test_restricting_the_attribute_names_converts_only_those(old_lut):
    path, _ = old_lut
    order_before = attribute_order(path)
    result = sa.repair_file(path, attribute_names={"spacing"}, min_age_seconds=600)
    assert result.status == "repaired" and len(result.conversions) == 5
    types = text_types(path)
    assert all(types[f"{axis}:spacing"] == "NC_STRING" for axis in sa.ORAC_AXIS_VARIABLES)
    assert types["optical_depth:units"] == "NC_CHAR" and types["optical_depth:long_name"] == "NC_CHAR"
    assert attribute_order(path) == order_before
    with Dataset(path) as ds:
        assert ds.variables["optical_depth"].getncattr("units") == "dimensionless"


def test_multiple_affected_attributes_across_variables_and_global_scope(old_lut):
    path, _ = old_lut
    audit = sa.audit_file(path)
    scopes = {c.variable for c in audit.conversions}
    assert None in scopes and len(scopes) >= 8
    assert sa.repair_file(path, min_age_seconds=600).status == "repaired"
    assert sa.audit_file(path).compliant


@pytest.mark.parametrize("kind", ["invalid_utf8", "user_defined_type"])
def test_unexpected_attribute_types_block_the_repair(old_lut, kind):
    path, _ = old_lut
    with nc.open_file(path, nc.NC_WRITE) as ncid:
        varid = nc.inq_varid(ncid, "optical_depth")
        if kind == "invalid_utf8":
            nc.put_att_text(ncid, varid, "odd", b"\xff\xfe not utf-8")
        else:
            typeid = nc.def_vlen(ncid, "vlen_int", nc.NC_INT)
            nc.put_att_vlen_int(ncid, varid, "odd", typeid, [1, 2, 3])
    os.utime(path, (OLD, OLD))
    digest = sha256(path)
    audit = sa.audit_file(path)
    assert audit.problems and "optical_depth:odd" in audit.problems[0]
    result = sa.repair_file(path, min_age_seconds=600)
    assert result.status == "failed" and "odd" in result.message
    assert sha256(path) == digest
    with pytest.raises(ValueError):
        sa.convert_attributes(path, audit)


def test_failed_conversion_leaves_the_original_untouched(old_lut, monkeypatch):
    path, _ = old_lut
    digest, mtime = sha256(path), path.stat().st_mtime

    def explode(*args, **kwargs):
        raise RuntimeError("simulated crash during nc_put_att_string")

    monkeypatch.setattr(sa, "convert_attributes", explode)
    result = sa.repair_file(path, min_age_seconds=600)
    assert result.status == "failed" and "original unchanged" in result.message
    assert sha256(path) == digest and path.stat().st_mtime == mtime
    assert not path.with_name(path.name + sa.TEMP_SUFFIX).exists()


def test_failed_verification_leaves_the_original_untouched(old_lut, monkeypatch):
    path, _ = old_lut
    digest = sha256(path)
    monkeypatch.setattr(sa, "verify_preserved",
                        lambda *a, **k: sa.VerificationReport(ok=False, differences=["simulated difference"]))
    result = sa.repair_file(path, min_age_seconds=600)
    assert result.status == "failed" and "simulated difference" in result.message
    assert sha256(path) == digest
    assert text_types(path)["optical_depth:spacing"] == "NC_CHAR"
    assert not path.with_name(path.name + sa.TEMP_SUFFIX).exists()


def test_verification_detects_tampered_data_and_metadata(old_lut, tmp_path):
    path, reference = old_lut
    audit = sa.audit_file(reference)
    tampered = tmp_path / "tampered.nc"
    shutil.copyfile(path, tampered)
    sa.convert_attributes(tampered, audit)
    assert sa.verify_preserved(reference, tampered, audit.conversions).ok
    with Dataset(tampered, "r+") as ds:
        v = ds.variables["T_dd"]
        v.set_auto_maskandscale(False)
        v[1, 0] = np.float32(np.nan)          # a NaN with a different payload
    report = sa.verify_preserved(reference, tampered, audit.conversions)
    assert not report.ok and any("T_dd" in d and "data bytes differ" in d for d in report.differences)
    shutil.copyfile(path, tampered)
    sa.convert_attributes(tampered, audit)
    with Dataset(tampered, "r+") as ds:
        ds.variables["optical_depth"].setncattr("valid_range", np.array([0.0, 999.0], dtype=np.float32))
    report = sa.verify_preserved(reference, tampered, audit.conversions)
    assert not report.ok and any("optical_depth:valid_range" in d for d in report.differences)
    shutil.copyfile(path, tampered)
    sa.convert_attributes(tampered, audit)
    with Dataset(tampered, "r+") as ds:
        ds.variables["optical_depth"].setncattr_string("spacing", "linear")
    report = sa.verify_preserved(reference, tampered, audit.conversions)
    assert not report.ok and any("optical_depth:spacing" in d and "conversion not as expected" in d
                                 for d in report.differences)


def test_recent_files_and_fresh_working_copies_are_skipped(old_lut):
    path, _ = old_lut
    digest = sha256(path)
    now = time.time()
    os.utime(path, (now, now))
    assert sa.repair_file(path, min_age_seconds=600).status == "skipped"
    os.utime(path, (OLD, OLD))
    temp = path.with_name(path.name + sa.TEMP_SUFFIX)
    temp.write_bytes(b"in progress")
    assert sa.repair_file(path, min_age_seconds=600).status == "skipped"
    assert temp.exists() and sha256(path) == digest
    os.utime(temp, (OLD, OLD))                # stale copy from an interrupted run
    result = sa.repair_file(path, min_age_seconds=600)
    assert result.status == "repaired" and not temp.exists()


def test_original_can_be_retained_as_a_hard_link(old_lut, tmp_path):
    path, reference = old_lut
    inode_before = path.stat().st_ino
    retention = tmp_path / "originals"
    result = sa.repair_file(path, min_age_seconds=600, original_dir=retention)
    assert result.status == "repaired"
    retained = retention / path.name
    assert result.original_retained_as == str(retained)
    assert retained.stat().st_ino == inode_before and path.stat().st_ino != inode_before
    assert sha256(retained) == sha256(reference)
    # A second repair must not overwrite a retained original.
    write_old_style_lut(path)
    assert sa.repair_file(path, min_age_seconds=600, original_dir=retention).status == "skipped"
    assert sha256(retained) == sha256(reference)


def test_slab_ranges_are_bounded_and_cover_the_leading_axis():
    ranges = list(sa.slab_ranges((100, 1000), 4, 10_000))
    assert ranges[0] == (0, 2) and ranges[-1] == (98, 100)
    assert all(stop - start <= 2 for start, stop in ranges)
    assert [start for start, _ in ranges] == list(range(0, 100, 2))
    assert list(sa.slab_ranges((3, 1000), 4, 100)) == [(0, 1), (1, 2), (2, 3)]   # one oversized row per slab
    assert list(sa.slab_ranges((), 4, 100)) == [(0, 1)]
    assert list(sa.slab_ranges((5,), 8, 16)) == [(0, 2), (2, 4), (4, 5)]


def test_verification_reads_large_variables_in_bounded_slabs(old_lut):
    path, reference = old_lut
    audit = sa.audit_file(reference)
    sa.convert_attributes(path, audit)
    report = sa.verify_preserved(reference, path, audit.conversions, slab_bytes=256)
    assert report.ok, report.differences
    assert report.slabs_compared > report.variables_compared
    assert sa.repair_file(path, min_age_seconds=600, slab_bytes=256).status in ("compliant", "repaired")


def test_command_line_dry_run_check_inventory_and_repair(old_lut, tmp_path, capsys):
    path, _ = old_lut
    digest = sha256(path)
    assert sa.main(["--dry-run", str(path)]) == 0
    assert sha256(path) == digest and "needs repair" in capsys.readouterr().out
    assert sa.main(["--check", "--quiet", str(path)]) == 3
    assert sha256(path) == digest
    inventory, logs = tmp_path / "inventory.tsv", tmp_path / "logs"
    assert sa.main(["--inventory", str(inventory), "--log-dir", str(logs), str(path)]) == 0
    rows = inventory.read_text().splitlines()
    assert rows[0].split("\t") == list(sa.INVENTORY_COLUMNS)
    fields = dict(zip(sa.INVENTORY_COLUMNS, rows[1].split("\t")))
    assert fields["status"] == "repaired" and fields["version"] == "25" and fields["required_conversions"] == "5"
    assert (logs / (path.name + ".repair.json")).exists()
    assert sa.main(["--check", str(path)]) == 0
    missing = tmp_path / "missing.nc"
    assert sa.main(["--dry-run", str(missing)]) == 1


def test_idl_written_v11_reference_is_compliant():
    if not V11_REFERENCE.exists():
        pytest.skip("authorised ORAC reference archive is unavailable")
    audit = sa.audit_file(V11_REFERENCE)
    assert audit.compliant
    assert audit.text_attribute_types == {"NC_STRING": 86}
    assert nc.read_string_attribute_like_orac(V11_REFERENCE, "optical_depth", "spacing") == "uneven_logarithmic"
