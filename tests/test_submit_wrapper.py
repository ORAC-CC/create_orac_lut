"""Shell-level tests of scripts/submit_oraclut_lut.sh with mocked hostname, ssh and sbatch.

No real ssh connection is opened and no SLURM job is submitted: ``hostname`` is a
fake on PATH, and the wrapper's ``ORACLUT_SSH`` / ``ORACLUT_SBATCH`` hooks point at
recording scripts.  The "dispatched" tests use a fake ssh that executes the remote
command locally under a fake atmlxint7 hostname, so the whole two-hop path runs.
"""

from __future__ import annotations

import json
import os
import shlex
import stat
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
WRAPPER = ROOT / "scripts" / "submit_oraclut_lut.sh"
SUBMIT = ROOT / "submit"
CONFIG = "configs/meteosat-10_seviri_cloud_liquid-water_stg_test.driver"


def _script(path: Path, body: str) -> Path:
    path.write_text("#!/bin/bash\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


@pytest.fixture
def mocks(tmp_path):
    """Return a factory building a mock environment for a given local hostname."""

    def build(local_host: str, *, ssh_mode: str = "execute", remote_host: str = "atmlxint7"):
        bindir = tmp_path / f"bin_{local_host}_{ssh_mode}"
        bindir.mkdir()
        record = tmp_path / f"record_{local_host}_{ssh_mode}"
        # hostname reads its answer from a file so the "remote" side can switch it.
        hostfile = tmp_path / f"hostname_{local_host}_{ssh_mode}"
        hostfile.write_text(local_host)
        _script(bindir / "hostname", f'cat "{hostfile}"\n')
        _script(bindir / "sbatch", (
            f'echo "call" >> "{record}.calls"\n'
            f'echo "SBATCH_ARGS" > "{record}.sbatch"\n'
            f'for a in "$@"; do printf "%s\\n" "$a" >> "{record}.sbatch"; done\n'
            f'echo "TMPDIR=${{TMPDIR-<unset>}} TMP=${{TMP-<unset>}} TEMP=${{TEMP-<unset>}}" >> "{record}.sbatch"\n'
            'echo "Submitted batch job 424242"\n'
        ))
        if ssh_mode == "execute":
            # Record the ssh argv, then behave like sshd: space-join the remote
            # arguments and hand the string to a shell on the "remote" host.
            _script(bindir / "ssh", (
                f'for a in "$@"; do printf "%s\\n" "$a" >> "{record}.ssh"; done\n'
                f'echo "{remote_host}" > "{hostfile}"\n'
                'while [[ "$1" == -o ]]; do shift 2; done\n'
                'shift  # host\n'
                f'export PATH="{bindir}:$PATH"\n'
                # The wrapper asks for a login shell so `module`/sbatch are on
                # PATH on the real atmlxint7; here a login profile would reset
                # PATH and hide the fake hostname, so run it as a plain shell.
                'cmd="$*"; exec bash -c "${cmd/bash -lc/bash -c}"\n'
            ))
        elif ssh_mode == "kerberos":
            _script(bindir / "ssh", (
                f'for a in "$@"; do printf "%s\\n" "$a" >> "{record}.ssh"; done\n'
                'echo "Permission denied (gssapi-keyex,gssapi-with-mic)." >&2\n'
                'exit 255\n'
            ))
        env = dict(os.environ)
        env["PATH"] = f"{bindir}:{env['PATH']}"
        env["ORACLUT_SBATCH"] = str(bindir / "sbatch")
        env["ORACLUT_SSH"] = str(bindir / "ssh")
        env["TMPDIR"] = "/tmp/user/27004"
        env["TMP"] = "/tmp/user/27004"
        env.pop("ORACLUT_DISPATCHED", None)
        return env, record

    return build


def _run(env, *args):
    return subprocess.run([str(WRAPPER), *args], cwd=ROOT, env=env, capture_output=True, text=True)


def _run_submit(env, *args, cwd=ROOT, command=SUBMIT):
    return subprocess.run([str(command), *args], cwd=cwd, env=env, capture_output=True, text=True)


def _sbatch_record(record: Path) -> list[str]:
    return (Path(str(record) + ".sbatch")).read_text().splitlines()


def _ssh_record(record: Path) -> list[str]:
    path = Path(str(record) + ".ssh")
    return path.read_text().splitlines() if path.exists() else []


def test_atmlxint7_calls_sbatch_directly_without_ssh(mocks):
    env, record = mocks("atmlxint7")
    result = _run(env, "--config", CONFIG, "--time=01:00:00")
    assert result.returncode == 0, result.stderr
    assert "Submitted batch job 424242" in result.stdout
    assert "Dispatching" not in result.stdout
    assert _ssh_record(record) == []
    lines = _sbatch_record(record)
    assert lines[1:4] == ["--time=01:00:00", "scripts/run_oraclut_lut.slurm", "--config"]
    assert lines[4] == CONFIG


def test_atmlxint5_dispatches_exactly_once_and_reaches_sbatch(mocks):
    env, record = mocks("atmlxint5")
    result = _run(env, "--config", CONFIG, "--time=12:00:00", "--mem=16G", "--cpus-per-task=1")
    assert result.returncode == 0, result.stderr
    assert "Dispatching ORAC LUT submission from atmlxint5 to atmlxint7" in result.stdout
    assert result.stdout.count("Dispatching") == 1
    assert "Submitting ORAC LUT job from atmlxint7" in result.stdout
    ssh = _ssh_record(record)
    assert ssh.count("atmlxint7") == 1, "exactly one ssh hop"
    assert "-o" in ssh and "BatchMode=yes" in ssh
    remote = ssh[-1]
    # ssh argv ends with the %q-quoted remote command handed to `bash -lc`.
    assert ssh[-3:-1] == ["bash", "-lc"]
    assert shlex.split(remote)[0].startswith(f"cd {ROOT} && ORACLUT_DISPATCHED=1 exec scripts/submit_oraclut_lut.sh")
    # Every resource and generator argument survives the hop and reaches sbatch.
    lines = _sbatch_record(record)
    assert lines[1:4] == ["--time=12:00:00", "--mem=16G", "--cpus-per-task=1"]
    assert lines[4] == "scripts/run_oraclut_lut.slurm"
    assert lines[5:7] == ["--config", CONFIG]
    # The sbatch boundary sees no inherited temporary-directory variables.
    assert lines[-1] == "TMPDIR=<unset> TMP=<unset> TEMP=<unset>"


def test_generator_arguments_and_spaces_survive_dispatch(mocks, tmp_path):
    env, record = mocks("atmlxint5")
    spaced = ROOT / "validation" / "generated" / "wrapper test dir"
    spaced.mkdir(parents=True, exist_ok=True)
    try:
        result = _run(
            env, "--time", "02:00:00", "--config", CONFIG, "--output", str(spaced),
            "--channels", "1, 9", "--dry-run", "--overwrite", "--job-name=my job",
        )
        assert result.returncode == 0, result.stderr
        lines = _sbatch_record(record)
        assert "--time" in lines and "02:00:00" in lines and "--job-name=my job" in lines
        assert str(spaced) in lines                     # single argv entry despite the space
        assert "1, 9" in lines and "--dry-run" in lines and "--overwrite" in lines
        # The remote command line quotes the spaced values rather than splitting them.
        remote = shlex.split(_ssh_record(record)[-1])[0]   # unquote the bash -lc argument
        assert shlex.split(remote.split("exec ", 1)[1])[1:] == [
            "--time", "02:00:00", "--config", CONFIG, "--output", str(spaced),
            "--channels", "1, 9", "--dry-run", "--overwrite", "--job-name=my job",
        ]
    finally:
        spaced.rmdir()


def test_missing_configuration_fails_before_any_dispatch(mocks):
    env, record = mocks("atmlxint5")
    result = _run(env, "--config", "configs/does_not_exist.driver")
    assert result.returncode == 2
    assert "configuration file not found" in result.stderr
    assert _ssh_record(record) == []
    assert not Path(str(record) + ".sbatch").exists()


def test_kerberos_failure_gives_kinit_guidance(mocks):
    env, record = mocks("atmlxint5", ssh_mode="kerberos")
    result = _run(env, "--config", CONFIG, "--time=12:00:00")
    assert result.returncode == 3
    assert "expired Kerberos ticket" in result.stderr
    assert "kinit" in result.stderr
    assert "Retry:  scripts/submit_oraclut_lut.sh --config" in result.stderr
    assert "--time=12:00:00" in result.stderr
    assert not Path(str(record) + ".sbatch").exists()


def test_dispatched_invocation_on_wrong_host_refuses_to_recurse(mocks):
    env, record = mocks("atmlxint5", ssh_mode="execute", remote_host="atmlxint6")
    result = _run(env, "--config", CONFIG)
    assert result.returncode != 0
    assert "landed on 'atmlxint6', not atmlxint7" in result.stderr
    assert _ssh_record(record).count("atmlxint7") == 1, "no second hop was attempted"
    assert not Path(str(record) + ".sbatch").exists()


def test_dispatched_marker_on_a_dispatch_host_is_refused(mocks):
    env, record = mocks("atmlxint5")
    env["ORACLUT_DISPATCHED"] = "1"
    result = _run(env, "--config", CONFIG)
    assert result.returncode == 2 and "refusing to dispatch again" in result.stderr
    assert _ssh_record(record) == []


@pytest.mark.parametrize("host", ["atmlxint4", "atmnode014", "somelaptop"])
def test_other_hosts_are_refused(mocks, host):
    env, record = mocks(host)
    result = _run(env, "--config", CONFIG)
    assert result.returncode == 2
    assert "only permitted from atmlxint7" in result.stderr
    assert _ssh_record(record) == [] and not Path(str(record) + ".sbatch").exists()


def test_help_and_scripts_parse():
    assert subprocess.run(["bash", "-n", str(WRAPPER)], capture_output=True).returncode == 0
    assert subprocess.run(["bash", "-n", str(ROOT / "scripts" / "run_oraclut_lut.slurm")], capture_output=True).returncode == 0
    result = subprocess.run([str(WRAPPER), "--help"], cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 0 and "kinit" in result.stdout and "atmlxint7" in result.stdout


RUN_FILE = "runs/meteosat-10_seviri_cloud_liquid-water_stg_test_ch01_ch09.run"


def test_run_file_is_dispatched_to_sbatch_unchanged(mocks):
    # IDL-structured path: the run file is the only generator argument and it
    # reaches the batch script exactly as given.
    env, record = mocks("atmlxint5")
    result = _run(env, RUN_FILE, "--time=01:30:00")
    assert result.returncode == 0, result.stderr
    lines = _sbatch_record(record)
    assert lines[1:4] == ["--time=01:30:00", "scripts/run_oraclut_lut.slurm", RUN_FILE]


def test_missing_run_file_fails_before_any_dispatch(mocks):
    env, record = mocks("atmlxint5")
    result = _run(env, "runs/does_not_exist.run")
    assert result.returncode == 2
    assert "run file not found" in result.stderr
    assert _ssh_record(record) == []
    assert not Path(str(record) + ".sbatch").exists()


def test_run_file_rejects_extra_generator_arguments(mocks):
    env, record = mocks("atmlxint7")
    result = _run(env, RUN_FILE, "--channels", "1")
    assert result.returncode == 2
    assert "takes no further generator arguments" in result.stderr
    assert not Path(str(record) + ".sbatch").exists()


def test_batch_script_routes_a_run_file_to_create_orac_luts():
    text = (ROOT / "scripts" / "run_oraclut_lut.slurm").read_text()
    assert '*.run)' in text
    assert '"$PYTHON" create_orac_luts.py "$1"' in text
    assert '"$PYTHON" -m oraclut.generate "$@"' in text
    assert '"$PYTHON" -m oraclut.master_run makerunfile --manifest-task' in text


def test_master_run_preflights_and_submits_one_array(mocks):
    env, record = mocks("atmlxint7")
    result = _run(env, "makerunfile", "--time=01:00:00")
    assert result.returncode == 0, result.stderr
    assert "Master expansion preflight" in result.stdout
    assert "liquid-water_old.mm" in result.stdout
    assert len((Path(str(record) + ".calls")).read_text().splitlines()) == 1
    lines = _sbatch_record(record)
    assert "--array=0-0" in lines
    assert "--master-manifest" in lines and "validation/slurm/manifests/" in lines[lines.index("--master-manifest") + 1]
    assert "--time=01:00:00" in lines


def test_submit_command_uses_fixed_master_and_default_resources(mocks, tmp_path):
    env, record = mocks("atmlxint7")
    result = _run_submit(env, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    lines = _sbatch_record(record)
    assert lines[1:4] == ["--time=12:00:00", "--mem=16G", "--cpus-per-task=1"]
    assert "--array=0-0" in lines
    manifest = Path(lines[lines.index("--master-manifest") + 1])
    assert str(manifest).startswith("validation/slurm/manifests/")
    assert json.loads(manifest.read_text())["master"] == "makerunfile"


def test_submit_command_allows_resource_overrides(mocks):
    env, record = mocks("atmlxint7")
    result = _run_submit(env, "--time=24:00:00", "--mem=32G")
    assert result.returncode == 0, result.stderr
    lines = _sbatch_record(record)
    assert "--time=24:00:00" in lines and "--mem=32G" in lines
    assert "--array=0-0" in lines


def test_submit_command_resolves_repository_through_symlink(mocks, tmp_path):
    env, record = mocks("atmlxint7")
    link_dir = tmp_path / "personal-bin"
    link_dir.mkdir()
    link = link_dir / "submit"
    link.symlink_to(SUBMIT)
    result = _run_submit(env, cwd=tmp_path, command=link)
    assert result.returncode == 0, result.stderr
    lines = _sbatch_record(record)
    assert lines[1:4] == ["--time=12:00:00", "--mem=16G", "--cpus-per-task=1"]
    assert lines[4] == "--array=0-0"
    assert lines[5] == "scripts/run_oraclut_lut.slurm"
    manifest = Path(lines[lines.index("--master-manifest") + 1])
    assert json.loads(manifest.read_text())["master"] == "makerunfile"
