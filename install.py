#!/usr/bin/env python3
"""Install (or remove) QRef in a Phenix installation."""
import argparse
import collections
import difflib
import glob
import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
import sys

BEGIN = '# QRef INSERT'
END = '# QRef END'
BACKUP_SUFFIX = '.qref-orig'
MANIFEST = '.qref-install.json'

BSD = 'cctbx, BSD 3-clause'
PROPRIETARY = 'Phenix, non-commercial source licence'

HERE = os.path.dirname(os.path.abspath(__file__))


class EditError(Exception):
    pass


# --- locating edit sites ---------------------------------------------------

# locate() pipes this to the cctbx.python of the installation being patched;
# the parsing happens there because the target files may use syntax this
# interpreter cannot read, and that interpreter can be Python 2.7
LOCATOR = r'''
from __future__ import print_function
import ast, sys


def fail(message):
    print('ERROR %s' % message)
    raise SystemExit(0)


def span(node):
    """First and last source line of a statement."""
    last = getattr(node, 'end_lineno', None)
    if last is None:                    # Python < 3.8 has no end_lineno
        last = max([getattr(n, 'lineno', node.lineno) for n in ast.walk(node)])
    return node.lineno, last


def emit(hits, want, what):
    """Print one SPAN per hit, after checking there are as many as expected."""
    if not hits:
        fail('no %s' % what)
    if want and len(hits) != want:
        fail('%d %s, expected %d' % (len(hits), what, want))
    for first, last in sorted(set(hits)):
        print('SPAN %d %d' % (first, last))


def body(cls, func):
    """Every node inside the named method."""
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == cls:
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name == func:
                    return list(ast.walk(item))
    fail('no method %s.%s' % (cls, func))


def called(node):
    """Trailing name of the function an Expr node calls: self.a.b() -> 'b'."""
    func = node.value.func
    return func.attr if isinstance(func, ast.Attribute) else getattr(func, 'id', '')


path, query, rest = sys.argv[1], sys.argv[2], sys.argv[3:]
try:
    tree = ast.parse(open(path).read())
except SyntaxError:
    fail('%s does not parse: %s' % (path, sys.exc_info()[1]))

if query == 'import':
    want = rest[0]
    for node in tree.body:
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == want and alias.asname is None:
                    raise SystemExit(0)         # already imported
    futures = [n for n in tree.body
               if isinstance(n, ast.ImportFrom) and n.module == '__future__']
    imports = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))]
    if futures:
        at = span(futures[-1])[1] + 1
    elif imports:
        at = imports[0].lineno
    elif tree.body:
        first = tree.body[0]
        is_docstring = (isinstance(first, ast.Expr)
                        and isinstance(getattr(first, 'value', None), ast.Str))
        at = span(first)[1] + 1 if is_docstring else first.lineno
    else:
        fail('%s is empty' % path)
    print('SPAN %d %d' % (at, at - 1))          # empty span: insert before `at`

elif query == 'stmt_call':
    cls, func = rest[0], rest[1]
    attr = rest[2] if len(rest) > 2 else None
    hits = [span(s) for s in body(cls, func)
            if isinstance(s, ast.Expr) and isinstance(s.value, ast.Call)
            and (attr is None or called(s) == attr)]
    emit(hits, 1 if attr else 0, 'call statements in %s.%s%s'
         % (cls, func, ' to ' + attr if attr else ''))

elif query == 'stmt_assign':
    cls, func, name = rest
    hits = [span(s) for s in body(cls, func) if isinstance(s, ast.Assign)
            and any(isinstance(t, ast.Name) and t.id == name for t in s.targets)]
    emit(hits, 0, 'assignments to %s in %s.%s' % (name, cls, func))

elif query == 'importfrom':
    part, name = rest
    hits = [span(n) for n in ast.walk(tree)
            if isinstance(n, ast.ImportFrom) and n.module and part in n.module
            and any(a.name == name for a in n.names)]
    emit(hits, 1, 'imports of %s from *%s*' % (name, part))

else:
    fail('unknown query %r' % query)
'''


def locate(layout, path, query):
    """Ask LOCATOR for the lines a query matches, as (first, last) pairs.

    It answers 'SPAN first last' per match, 1-based and inclusive.  No pairs
    means there is nothing to do, and 'ERROR message' at exit 0 means LOCATOR
    could not answer; a non-zero exit is cctbx.python itself failing.
    """
    if layout.python is None:
        raise EditError(f'no cctbx.python for this installation, so {path} '
                        f'cannot be parsed')
    argv = [layout.python, '-', path] + [str(a) for a in query]
    proc = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    out = proc.communicate(LOCATOR.encode('utf-8'))[0]
    if not isinstance(out, str):
        out = out.decode('utf-8', 'replace')
    if proc.returncode != 0:
        raise EditError(f'parsing {path} failed:\n{out.strip()}')
    spans = []
    for line in out.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == 'ERROR':
            raise EditError(line.strip()[len('ERROR '):])
        spans.append((int(fields[1]), int(fields[2])))
    return spans


def indent_of(lines, number):
    line = lines[number - 1] if 0 < number <= len(lines) else ''
    return line[:len(line) - len(line.lstrip())] if line.strip() else ''


# an edit is (query, block before each match, block after it); either can be None
def ensure_import(module):
    return (('import', module), [f'import {module} # QRef'], None)


def lock(name):
    """Create a lock file around the wrapped statement; no finally, so an exception leaves it."""
    return ([BEGIN, f"with open('{name}', 'w'):", '  pass'],
            [f"os.remove('{name}')", END])


# a planned insertion: the 1-based line it goes before, its rendered source, and
# whether it opens a statement rather than closing the one before it
Block = collections.namedtuple('Block', 'at lines opening')


def rendered(block, indent):
    """A block of source, indented to sit beside the statement it brackets."""
    return [(indent + text if text.strip() else '') + '\n' for text in block]


def plan(layout, path, lines, edits):
    """Blocks to insert, before each match's first line and after its last."""
    blocks = []
    for query, before, after in edits:
        for first, last in locate(layout, path, query):
            indent = indent_of(lines, first)
            if before:
                blocks.append(Block(first, rendered(before, indent), True))
            if after:
                blocks.append(Block(last + 1, rendered(after, indent), False))
    return blocks


def splice(lines, blocks):
    """Insert the planned blocks, working up the file so the line numbers hold.

    Two blocks can want the same line, where one statement's closing block meets
    the next statement's opening block.  The opening block goes in first, so it
    ends up second.
    """
    out = list(lines)
    for block in sorted(blocks, key=lambda b: (-b.at, not b.opening)):
        out[block.at - 1:block.at - 1] = block.lines
    return out


# --- what to insert --------------------------------------------------------

def qref_hook(guard):
    """The call to QRef; the unused guard is emitted commented out."""
    lock_test = ("if not os.path.exists('qm.lock') and "
                 "(os.path.exists('xyz_reciprocal.lock') or os.path.exists('xyz.lock')):")
    plain_test = "if not os.path.exists('qm.lock'):"
    active, idle = (lock_test, plain_test) if guard == 'locks' else (plain_test, lock_test)
    return [
        BEGIN,
        "if compute_gradients is True and os.path.exists('qref.dat'):",
        '  # ' + idle,
        '  ' + active,
        '    from qref import qref',
        '    self.gradients, self.residual_sum = qref.run(sites_cart=sites_cart, '
        'mm_gradients=self.gradients, mm_residual_sum=self.residual_sum)',
        END,
    ]


# dumps the pdb_interpretation_params argument; QRef reads settings.pickle back
SETTINGS_PICKLE = [
    BEGIN,
    "if not os.path.exists('settings.pickle'):",
    "  with open('settings.pickle', 'wb') as file:",
    '    pickle.dump(pdb_interpretation_params, file)',
    END,
]

# real_space_refine leaves the output model translated by -shift_cart when
# sort_atoms = False; this prints the translation that undoes it, at the end
# of refinement in the .log file
# anchored on the shift call in shift_model_back, past the assert that sets
# shift_cart()
SHIFT_HINT = [
    BEGIN,
    '# real_space_refine does not shift the model back to the original coordinate',
    '# frame when sort_atoms = False; apply this translation to the output model',
    '# if it comes out displaced.',
    'print()',
    "print('QRef: if the output model is not shifted back to the original frame, apply')",
    'print(\'QRef:   phenix.pdbtools translate="\' + '
    "' '.join([str(round(x, 4)) for x in self.shift_cart()]) + '\"')",
    'print()',
    END,
]

class Target(collections.namedtuple('Target', 'root rel edits licence')):
    """One file QRef edits: where it lives, what goes into it, whose it is."""

    @property
    def proprietary(self):
        return self.licence is PROPRIETARY

    def mode_for(self, args):
        """What to do with this file: skip, report, diff or apply."""
        if self.proprietary and args.phenix_edits != 'apply':
            return args.phenix_edits            # 'skip' or 'report'
        return 'diff' if args.dry_run else 'apply'

    def report(self, report):
        """Print where the edits go, for someone applying them by hand."""
        path = report.layout.path(self)
        lines = read_lines(path)
        if any(BEGIN in line for line in lines):
            report.entry(self.rel, 'already patched, nothing to do')
            return
        detail = ['file ' + path]
        for block in sorted(plan(report.layout, path, lines, self.edits),
                            key=lambda b: -b.at):
            detail.append(f'insert before line {block.at}:')
            detail += ['    ' + line.rstrip('\n') for line in block.lines]
        report.entry(self.rel, 'edits below, apply by hand', detail)

    def patch(self, report, force):
        """Insert QRef's blocks.  True if the file would be written."""
        layout = report.layout
        path = layout.path(self)
        if not os.path.exists(path):
            raise EditError(f'{self.rel}: missing ({path})')
        installed = read_lines(path)
        patched = any(BEGIN in line for line in installed)
        backup = path + BACKUP_SUFFIX

        if patched and not force:
            report.entry(self.rel, 'already patched, skipped')
            return False

        original, source, refreshed = installed, path, None
        if patched:
            # start from the backup, so that --force does not nest insertions
            if not os.path.exists(backup):
                raise EditError(f'{self.rel}: patched, but {backup} is gone; '
                                f'reinstall Phenix or revert the file by hand')
            original, source = read_lines(backup), backup
        elif not report.dry_run:
            refreshed = check_backup(path, backup)

        blocks = plan(layout, source, original, self.edits)
        lines = splice(original, blocks)
        status = (report.verb('would patch', 'patched') + ', '
                  + plural(len(blocks), 'insertion'))

        if report.dry_run:
            # for the two Phenix/refinement files, diff without context: only
            # what QRef adds
            detail = []
            if self.proprietary:
                detail.append(f'diff shown without context: {self.licence}')
            diff = difflib.unified_diff(original, lines, 'a/' + self.rel,
                                        'b/' + self.rel,
                                        n=0 if self.proprietary else 2)
            detail += [line.rstrip('\n') for line in diff]
            report.entry(self.rel, status, detail)
            return True

        commit(layout, path, lines, backup)
        report.entry(self.rel, status, [refreshed] if refreshed else [])
        return True

    def restore(self, report):
        """Put the backup back.  True if the file would be restored."""
        path = report.layout.path(self)
        backup = path + BACKUP_SUFFIX
        if not os.path.exists(backup):
            report.entry(self.rel, 'no backup, left alone')
            return False
        if report.dry_run:
            report.entry(self.rel, 'would restore')
            return True
        shutil.copy2(backup, path)
        os.remove(backup)
        drop_pycache(path)
        report.entry(self.rel, 'restored')
        return True


def patches(guard, shift_hint=False):
    """The files QRef edits, keyed by the name used in the output."""
    model_edits = [ensure_import('os'), ensure_import('pickle'),
                   (('importfrom', 'geometry_restraints', 'quantum_interface'),
                    SETTINGS_PICKLE, None)]
    if shift_hint:
        model_edits.append((('stmt_call', 'manager', 'shift_model_back',
                             'shift_model_and_set_crystal_symmetry'),
                            SHIFT_HINT, None))
    spec = {
        'energies.py': Target('cctbx', 'cctbx/geometry_restraints/energies.py', [
            ensure_import('os'),
            (('stmt_call', 'energies', '__init__',
              'finalize_target_and_gradients'), qref_hook(guard), None),
        ], BSD),
        'model.py': Target('cctbx', 'mmtbx/model/model.py', model_edits, BSD),
    }
    if guard != 'locks':
        return spec
    before, after = lock('xyz_reciprocal.lock')
    # minimized is assigned in two branches; each gets its own lock
    spec['xyz_reciprocal_space.py'] = Target(
        'phenix', 'phenix/refinement/xyz_reciprocal_space.py', [
            ensure_import('os'),
            (('stmt_assign', 'run_all', 'run_lbfgs', 'minimized'), before, after),
        ], PROPRIETARY)
    before, after = lock('xyz.lock')
    spec['macro_cycle_real_space.py'] = Target(
        'phenix', 'phenix/refinement/macro_cycle_real_space.py', [
            ensure_import('os'),
            (('stmt_call', 'run', 'refine_xyz'), before, after),
        ], PROPRIETARY)
    return spec


# --- the installation being written to -------------------------------------

class Layout(object):
    """Where a Phenix installation keeps its Python packages and interpreter."""

    def __init__(self, prefix):
        self.prefix = os.path.abspath(prefix)
        modules = os.path.join(self.prefix, 'modules')
        if os.path.isdir(os.path.join(modules, 'cctbx_project')):
            self.kind = 'modules'
            self.roots = {'cctbx': os.path.join(modules, 'cctbx_project'),
                          'phenix': os.path.join(modules, 'phenix')}
        else:
            self.kind = 'site-packages'
            root = self._site_packages()
            self.roots = {'cctbx': root, 'phenix': root}
        self.package_parent = self.roots['cctbx'] if self.kind == 'site-packages' else modules
        self.python = self._python()

    def _site_packages(self):
        # realpath() collapses the pythonX.Y -> pythonX.YZ symlink some installs ship
        found = sorted(set(
            os.path.realpath(path)
            for path in glob.glob(os.path.join(self.prefix, 'lib', 'python*',
                                               'site-packages'))
            if os.path.isdir(os.path.join(os.path.realpath(path), 'cctbx'))))
        if not found:
            raise EditError(f'no site-packages containing cctbx under '
                            f'{self.prefix} -- is this a Phenix installation?')
        if len(found) > 1:
            listed = '\n  '.join(found)
            raise EditError(f'several candidate site-packages under '
                            f'{self.prefix}:\n  {listed}')
        return found[0]

    def _python(self):
        for rel in ('phenix_bin/cctbx.python', 'build/bin/cctbx.python'):
            candidate = os.path.join(self.prefix, rel)
            if os.path.exists(candidate):
                return candidate
        return None

    def version(self):
        for name in ('phenix_env.sh', 'build/phenix_env.sh'):
            path = os.path.join(self.prefix, name)
            if os.path.exists(path):
                with open(path) as handle:
                    found = re.search(r'PHENIX_VERSION=(\S+)', handle.read())
                if found:
                    return found.group(1)
        return 'unknown'

    def path(self, target):
        return os.path.join(self.roots[target.root], target.rel)


# --- listing ---------------------------------------------------------------

COLUMN = 45                 # status column; the widest address is 43 characters


def under(path, base):
    """`path` relative to `base`, or None if it is not inside it.

    Both are tried as written and as resolved, so a symlink on either side does
    not hide the relationship.
    """
    for full in (os.path.abspath(path), os.path.realpath(path)):
        for root in (os.path.abspath(base), os.path.realpath(base)):
            root = root.rstrip(os.sep) + os.sep
            if full.startswith(root):
                return full[len(root):].replace(os.sep, '/')
    return None


def plural(count, word):
    return f'{count} {word}' + ('' if count == 1 else 's')


class Report(object):
    """The listing on screen, and whether this run writes anything.

    Carried by everything that reports what it did, so that the installation
    being written to, the dry-run flag and the output all travel together.
    """

    def __init__(self, layout, dry_run):
        self.layout = layout
        self.dry_run = dry_run
        self.blank = True               # to suppress a blank after a blank

    def line(self, message=''):
        if not message and self.blank:
            return
        print(message)
        self.blank = not message

    def heading(self, text):
        self.line()
        self.line(text)

    def field(self, name, value):
        self.line(f'{name:<11}{value}')

    def entry(self, address, status='', detail=()):
        """One listing line, the address and its status, with details under it.

        An entry whose detail runs to several lines is followed by a blank line,
        so that neighbouring blocks do not run together.
        """
        if not status:
            self.line('  ' + address)
        elif len(address) < COLUMN:
            self.line(f'  {address:<{COLUMN}}{status}')
        else:
            self.line('  ' + address)
            self.line('  ' + ' ' * COLUMN + status)
        detail = list(detail)
        while detail and not detail[-1].strip():
            detail.pop()
        for text in detail:
            self.line(('    ' + text).rstrip())
        if len(detail) > 1:
            self.line()

    def path(self, path):
        """A path as the listing addresses it: relative to its base, else absolute."""
        layout = self.layout
        for base in (layout.roots['cctbx'], layout.roots['phenix'], layout.prefix):
            relative = under(path, base)
            if relative is not None:
                return relative
        return os.path.abspath(path)

    def verb(self, future, past):
        """The word for something this run would do, or has done."""
        return future if self.dry_run else past

    def outcome(self):
        return 'dry run, nothing written' if self.dry_run else 'done'


# --- reading and writing ---------------------------------------------------

def read_lines(path):
    with open(path) as handle:
        return handle.readlines()


def write_lines(path, lines):
    with open(path, 'w') as handle:
        handle.writelines(lines)


def digest(path):
    with open(path, 'rb') as handle:
        return hashlib.md5(handle.read()).hexdigest()


def verify(layout, candidate, shown_as):
    """Byte-compile with the cctbx.python of the installation being patched."""
    proc = subprocess.Popen([layout.python, '-m', 'py_compile', candidate],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out = proc.communicate()[0]
    if proc.returncode != 0:
        if not isinstance(out, str):
            out = out.decode('utf-8', 'replace')
        out = out.replace(candidate, shown_as)
        raise EditError(f'the patched form of {shown_as} does not compile:\n'
                        f'{out}\nnothing was written; the installed file is '
                        f'unchanged')


def drop_pycache(path):
    """Remove any stale .pyc for path, in both Python 3 and 2 layouts."""
    directory, base = os.path.split(path)
    stem = os.path.splitext(base)[0]
    for stale in glob.glob(os.path.join(directory, '__pycache__', stem + '.*.pyc')):
        os.remove(stale)
    legacy = os.path.join(directory, stem + '.pyc')
    if os.path.exists(legacy):
        os.remove(legacy)


def commit(layout, path, lines, backup):
    """Write lines to path, via a temp file that has to byte-compile first."""
    directory, base = os.path.split(path)
    temp = os.path.join(directory, os.path.splitext(base)[0] + '.qref-new.py')
    write_lines(temp, lines)
    try:
        verify(layout, temp, path)
    except EditError:
        drop_pycache(temp)
        os.remove(temp)
        raise
    drop_pycache(temp)
    if not os.path.exists(backup):
        shutil.copy2(path, backup)
    mode = stat.S_IMODE(os.stat(path).st_mode)
    getattr(os, 'replace', os.rename)(temp, path)
    os.chmod(path, mode)
    drop_pycache(path)


def check_backup(path, backup):
    """Refresh a backup that no longer matches an unpatched file."""
    if os.path.exists(backup) and read_lines(backup) != read_lines(path):
        shutil.copy2(path, backup)
        return 'backup refreshed, the installed file had changed'
    return None


# --- patching --------------------------------------------------------------

# --- the package, the scripts and the templates ----------------------------

def install_package(report, checkout, link):
    """Put the qref package where this installation will import it from."""
    layout = report.layout
    source = os.path.join(checkout, 'qref')
    if not os.path.isdir(source):
        raise EditError(f'no qref package at {source}; pass --source <checkout>')
    package = os.path.join(layout.package_parent, 'qref')
    pth = os.path.join(layout.package_parent, 'qref.pth')
    if link and layout.kind != 'site-packages':
        raise EditError('--link needs the site-packages layout; the modules '
                        'layout has no .pth support')
    if link:
        destination = pth
        status = report.verb('would link', 'linked')
        detail = ['points at ' + checkout]
    else:
        destination = package
        status = report.verb('would install', 'installed')
        detail = ['from ' + source]
    report.entry(report.path(destination), status, detail)
    if report.dry_run:
        return destination
    if os.path.isdir(package):
        shutil.rmtree(package)
    if os.path.exists(pth):
        os.remove(pth)
    if link:
        write_lines(pth, [checkout + '\n'])
    else:
        shutil.copytree(source, package,
                        ignore=shutil.ignore_patterns('__pycache__', '*.pyc', '*.pyo'))
    return destination

class Extra(collections.namedtuple('Extra', 'option shebang plain')):
    """Files installed outside the package: the scripts and the templates.

    `option` is the command line option naming where they go, `shebang` the
    files whose first line is rewritten, `plain` those copied as they are.
    """

    def install(self, report, checkout, dest, pin=False):
        """Copy the files into `dest`.  Returns the paths written.

        With `pin`, the shebang of the files in `shebang` is set to this
        installation's cctbx.python; otherwise they keep `#!/usr/bin/env
        cctbx.python` and follow whichever Phenix is sourced when they run.

        utils.py travels with the scripts: they all import it, and python only
        puts a script's own directory on sys.path.
        """
        layout = report.layout
        source = os.path.join(checkout, self.option)
        if pin and self.shebang and layout.python is None:
            raise EditError(f'no cctbx.python under {layout.prefix}, cannot pin '
                            f'shebangs')
        if os.path.exists(dest) and not os.path.isdir(dest):
            raise EditError(f'{dest} exists and is not a directory')
        if not report.dry_run and not os.path.isdir(dest):
            os.makedirs(dest)
        written = []
        for name in self.shebang + self.plain:
            src, out = os.path.join(source, name), os.path.join(dest, name)
            if not os.path.exists(src):
                report.entry(name, f'missing in {source}, skipped')
                continue
            lines, detail = None, []
            if pin and name in self.shebang:
                lines = read_lines(src)
                had_shebang = lines and lines[0].startswith('#!')
                old = lines[0].rstrip('\n') if had_shebang else '<none>'
                detail = ['was ' + old]
                body = lines[1:] if had_shebang else lines
                lines = ['#!' + layout.python + '\n'] + body
                status = report.verb('would pin shebang', 'shebang pinned')
            else:
                status = report.verb('would copy', 'copied')
            report.entry(name, status, detail)
            if not report.dry_run:
                if lines is None:
                    shutil.copy2(src, out)
                else:
                    write_lines(out, lines)
                if name in self.shebang:
                    os.chmod(out, 0o755)    # whatever the checkout has
            written.append(out)
        return written


EXTRAS = (Extra('scripts', ('qref_prep.py', 'prep_geo_opt_qm_constrained.py',
                            'sort_pdb.py', 'change_occ_pdb.py'), ('utils.py',)),
          Extra('templates', (), ('junctfactor', 'qm_1.inp')))


# --- the manifest ----------------------------------------------------------

def manifest_path(layout):
    return os.path.join(layout.prefix, MANIFEST)


def write_manifest(report, files, directories):
    """Record what was installed outside the python package directories.

    --scripts and --templates write into directories the caller chose, which
    --uninstall has no other way to find.  Records are merged with any already
    there, so installing scripts and templates in separate runs still leaves
    both removable.  The digest lets --uninstall skip anything edited since.
    """
    if report.dry_run or not (files or directories):
        return
    path = manifest_path(report.layout)
    known = {}
    directories = set(directories)
    if os.path.exists(path):
        with open(path) as handle:
            old = json.load(handle)
        known = dict((item['path'], item) for item in old.get('files', []))
        directories |= set(old.get('directories', []))
    for name in files:
        known[name] = {'path': name, 'md5': digest(name)}
    with open(path, 'w') as handle:
        json.dump({'files': sorted(known.values(), key=lambda e: e['path']),
                   'directories': sorted(directories)},
                  handle, indent=2, sort_keys=True)
    report.heading('Manifest')
    report.entry(report.path(path), 'written, read back by --uninstall')


def remove_recorded(report, records, directory):
    """Remove the recorded files in one directory, and the directory if empty."""
    for record in sorted(records, key=lambda r: r['path']):
        name = os.path.basename(record['path'])
        if not os.path.exists(record['path']):
            report.entry(name, 'gone already')
        elif digest(record['path']) != record.get('md5'):
            report.entry(name, 'edited since install, left alone')
        elif report.dry_run:
            report.entry(name, 'would remove')
        else:
            os.remove(record['path'])
            drop_pycache(record['path'])    # left behind by running the script
            report.entry(name, 'removed')
    if not os.path.isdir(directory):
        return
    cache = os.path.join(directory, '__pycache__')
    if os.path.isdir(cache) and not os.listdir(cache) and not report.dry_run:
        os.rmdir(cache)
    if os.listdir(directory):
        return
    if report.dry_run:
        report.line('    directory would be left empty')
    else:
        os.rmdir(directory)
        report.line('    directory left empty, removed')


def remove_installed(report):
    """Undo what write_manifest recorded."""
    path = manifest_path(report.layout)
    if not os.path.exists(path):
        report.entry(report.path(path), 'not found',
                     ['scripts and templates left alone'])
        return
    with open(path) as handle:
        recorded = json.load(handle)
    report.entry(report.path(path),
                 report.verb('would be read', 'read, then removed'))
    grouped = {}
    for record in recorded.get('files', []):
        grouped.setdefault(os.path.dirname(record['path']), []).append(record)
    for directory in sorted(grouped):
        report.heading(f'Installed in {directory}')
        removable = directory in recorded.get('directories', [])
        remove_recorded(report, grouped[directory],
                        directory if removable else '')
    if not report.dry_run:
        os.remove(path)


# --- command line ----------------------------------------------------------

def parse_args(argv):
    parser = argparse.ArgumentParser(
        description='Install QRef into a Phenix installation.')
    parser.add_argument('--phenix', default=os.environ.get('PHENIX'),
                        help='Phenix prefix (default: $PHENIX)')
    parser.add_argument('--guard', choices=['locks', 'always'], default='locks',
                        help='"locks" runs QRef only while coordinates are '
                             'refined, which needs the two Phenix/refinement '
                             'files to set the lock (default); "always" runs it '
                             'whenever qref.dat is present and qm.lock is not, '
                             'leaving those two files untouched.')
    parser.add_argument('--phenix-edits', choices=['apply', 'report', 'skip'],
                        default='apply',
                        help='What to do with the two Phenix/refinement files: '
                             'edit in place (default), print the located edits '
                             'so you can apply them yourself, or leave alone.')
    parser.add_argument('--shift-hint', action='store_true',
                        help='Also print the translation needed to correct the '
                             'output of phenix.real_space_refine, which does not '
                             'shift the model back to the original coordinate '
                             'frame when sort_atoms = False.')
    parser.add_argument('--scripts', metavar='DIR',
                        help='Also install scripts/ here, shebangs rewritten to '
                             'this installation (put DIR on your PATH).')
    parser.add_argument('--pin-shebang', action='store_true',
                        help='if specified sets the shebang of the installed '
                             'scripts to the cctbx.python of this installation, '
                             'instead of leaving them to follow whichever '
                             'Phenix is sourced when they run')
    parser.add_argument('--templates', metavar='DIR',
                        help='Also copy templates/ here.')
    parser.add_argument('--link', action='store_true',
                        help='Point a .pth at the checkout instead of copying '
                             'the QRef package, so that several installations '
                             'share one source.')
    parser.add_argument('--source', default=HERE,
                        help='QRef checkout (default: the one this script is in).')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would change and touch nothing.')
    parser.add_argument('--force', action='store_true',
                        help='Re-apply over an existing install, starting from '
                             f'the {BACKUP_SUFFIX} backups.')
    parser.add_argument('--uninstall', action='store_true',
                        help='Restore the backups and remove what was installed.')
    return parser.parse_args(argv)


class Installer(object):
    """One run: what it is doing, where, and what it has done so far."""

    def __init__(self, report, args):
        self.report = report
        self.args = args
        self.layout = report.layout
        self.patched = []               # files written, named for the error note

    def header(self):
        report, layout, args = self.report, self.layout, self.args
        title = 'QRef uninstall' if args.uninstall else 'QRef install'
        if report.dry_run:
            title += ', dry run -- nothing is written'
        report.line(title)
        report.field('prefix', f'{layout.prefix}  (Phenix {layout.version()})')
        report.field('layout', layout.kind)
        roots = sorted(set(layout.roots.values()))
        report.field('packages', ', '.join(report.path(root) for root in roots))
        report.field('python', report.path(layout.python) if layout.python
                     else '<not found>')
        if not args.uninstall:
            report.field('source', args.source)
            report.field('guard', args.guard)
        report.line()

    def guard_note(self):
        self.report.line('note: --guard locks relies on the lock files that only '
                         'the two Phenix/refinement')
        self.report.line('      edits create, so QRef stays inactive until those '
                         'edits are in place.')
        self.report.line()

    def targets(self):
        """The files to work through, in the order the listing shows them."""
        spec = patches(self.args.guard, self.args.shift_hint)
        return sorted(spec.values(), key=lambda t: t.rel)

    def install(self):
        report, args = self.report, self.args
        if args.guard == 'locks' and args.phenix_edits != 'apply':
            self.guard_note()

        report.line('Module')
        added = [install_package(report, args.source, args.link)]

        heading = 'Patches'
        if not report.dry_run:
            heading += f' (originals kept as <file>{BACKUP_SUFFIX})'
        report.heading(heading)
        by_hand = 0
        for target in self.targets():
            mode = target.mode_for(args)
            if mode == 'skip':
                report.entry(target.rel, f'left alone ({target.licence})')
            elif mode == 'report':
                target.report(report)
                by_hand += 1
            elif target.patch(report, args.force):
                self.patched.append(target.rel)

        files, directories = [], []
        for extra in EXTRAS:
            dest = getattr(args, extra.option)
            if not dest:
                continue
            report.heading(f'{extra.option.capitalize()} -> {dest}')
            files += extra.install(report, args.source, dest,
                                   pin=args.pin_shebang)
            directories.append(os.path.abspath(dest))
        write_manifest(report, files, directories)

        report.line()
        would = 'would be ' if report.dry_run else ''
        parts = [f"{plural(len(self.patched), 'file')} {would}patched",
                 f'{len(added) + len(files)} {would}installed']
        if by_hand:
            parts.append(f"{plural(by_hand, 'file')} left for you")
        report.line(f"{report.outcome()}: {', '.join(parts)}")

    def uninstall(self):
        """Restore every backup and remove everything an install put in place."""
        report = self.report
        restored = 0
        report.line('Patches')
        # always the full set, whatever --guard and --phenix-edits chose at install
        for target in sorted(patches('locks', shift_hint=True).values(),
                             key=lambda t: t.rel):
            if target.restore(report):
                restored += 1
        report.heading('Module')
        for name in ('qref', 'qref.pth'):
            leftover = os.path.join(self.layout.package_parent, name)
            if not os.path.exists(leftover):
                continue
            if report.dry_run:
                report.entry(report.path(leftover), 'would remove')
            else:
                if os.path.isdir(leftover):
                    shutil.rmtree(leftover)
                else:
                    os.remove(leftover)
                report.entry(report.path(leftover), 'removed')
        report.heading('Manifest')
        remove_installed(report)
        report.line()
        report.line(f"{report.outcome()}: {plural(restored, 'file')} "
                    f"{report.verb('would be restored', 'restored')}")

    def fail(self, message, code):
        sys.stdout.flush()
        print(f'error: {message}', file=sys.stderr)
        if self.patched:
            print(f"note: {', '.join(self.patched)} already patched; run "
                  f"--uninstall to back that out", file=sys.stderr)
        return code


def main(argv=None):
    args = parse_args(argv)
    if not args.phenix:
        print('error: no Phenix prefix; pass --phenix or source phenix_env.sh',
              file=sys.stderr)
        return 2
    try:
        layout = Layout(args.phenix)
    except EditError as error:
        print(f'error: {error}', file=sys.stderr)
        return 2

    installer = Installer(Report(layout, args.dry_run), args)
    installer.header()
    try:
        if args.uninstall:
            installer.uninstall()
        else:
            installer.install()
    except (EditError, IOError, OSError) as error:
        return installer.fail(error, 1)
    return 0


if __name__ == '__main__':
    sys.exit(main())
