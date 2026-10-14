"""Tests for install.py that write to a temporary directory only."""
import os
import subprocess
import sys
from types import SimpleNamespace

import install

EXTRAS = dict((extra.option, extra) for extra in install.EXTRAS)


def report(tmp_path):
    """A Report for an installation whose cctbx.python is the one running the
    tests, so that pinned scripts can be run."""
    prefix = str(tmp_path)
    layout = SimpleNamespace(prefix=prefix, python=sys.executable,
                             roots={'cctbx': prefix, 'phenix': prefix})
    return install.Report(layout, dry_run=False)


def test_linked_scripts_run_from_the_checkout(tmp_path):
    """With --link, the scripts are wrappers that run the scripts in the
    checkout, so that pulling the checkout updates them; utils.py is left out,
    as they find it beside them in the checkout."""
    scripts = EXTRAS['scripts']
    dest = str(tmp_path / 'bin')
    scripts.install(report(tmp_path), install.HERE, dest, pin=True, link=True)

    assert sorted(os.listdir(dest)) == sorted(scripts.shebang)
    wrapper = os.path.join(dest, 'qref_prep.py')
    with open(wrapper) as handle:
        assert os.path.join(install.HERE, 'scripts', 'qref_prep.py') in handle.read()
    run = subprocess.run([wrapper, '--version'], stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, universal_newlines=True)
    assert run.returncode == 0, run.stdout
    assert run.stdout.splitlines()[-1].startswith('qref_prep.py ')


def test_uninstall_removes_linked_scripts(tmp_path):
    """The manifest records the wrappers, so --uninstall removes them along with
    the templates, which --link still copies."""
    out = report(tmp_path)
    dest = str(tmp_path / 'bin')
    files = []
    for extra in install.EXTRAS:
        files += extra.install(out, install.HERE, dest, pin=True, link=True)
    with open(os.path.join(dest, 'junctfactor')) as copied, \
            open(os.path.join(install.HERE, 'templates', 'junctfactor')) as source:
        assert copied.read() == source.read()

    install.write_manifest(out, files, [dest])
    install.remove_installed(out)
    assert not os.path.exists(dest)
