"""Smoke test: import + version + CLI invocation."""

from __future__ import annotations

from typer.testing import CliRunner


def test_import_and_version() -> None:
    import pd_target_credentialing

    assert pd_target_credentialing.__version__ == "0.1.0"


def test_cli_version() -> None:
    from pd_target_credentialing.cli import app

    runner = CliRunner()
    result = runner.invoke(app, ["version"])
    assert result.exit_code == 0
    assert "pd-target-credentialing" in result.stdout.lower() or "0.1.0" in result.stdout
