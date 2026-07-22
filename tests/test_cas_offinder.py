from pathlib import Path

from otp import cas_offinder


def test_ensure_cas_offinder_uses_cmake_policy_compatibility_flag(tmp_path, monkeypatch):
    repo_dir = tmp_path / "cas-offinder"
    repo_dir.mkdir()
    monkeypatch.setattr(cas_offinder, "DATA_DIR", tmp_path)
    monkeypatch.setattr(cas_offinder.shutil, "which", lambda _: None)

    original_path = cas_offinder.Path
    monkeypatch.setattr(
        cas_offinder,
        "Path",
        lambda value: tmp_path / "docker-bin" if value == "/usr/local/bin/cas-offinder" else original_path(value),
    )

    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        if command == ["make"]:
            (Path(kwargs["cwd"]) / "cas-offinder").touch()

    monkeypatch.setattr(cas_offinder.subprocess, "run", fake_run)

    assert cas_offinder.ensure_cas_offinder() == str(repo_dir / "build" / "cas-offinder")
    assert calls[0][0] == [
        "cmake",
        "-G",
        "Unix Makefiles",
        "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
        "..",
    ]
