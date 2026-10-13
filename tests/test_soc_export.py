"""Regression coverage for the SOC step of a properties/data run.

The spin-orbit coupling that reaches a dynamics driver leaves the engine twice:
as the ``socme`` table of singlet-triplet couplings, and as one set of energies
and gradients per spin manifold.  Both are checked against the same run's own
log, because the native kernel prints the very table the file must reproduce --
the transpose, the alpha^2/2 conversion and the row order are all Python-side
choices with no other witness.
"""

import re
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
SOC_DECK = EXAMPLES / "SOC" / "H2O_BHHLYP_MRSF_PROP_SOC.inp"
SOC_DECK_STANDALONE = EXAMPLES / "SOC" / "H2O_BHHLYP_SOC.inp"
PLAIN_PROP_DECK = EXAMPLES / "other" / "h2o_rohf_mrsf-s_6-31g_prop.inp"

#: One row of the native singlet-triplet table: labels, then Re(1e) Im(1e)
#: Re(2e) Im(2e) Absolute in cm-1.  The Ms=0 label carries an internal space
#: ("(Ms= 0)"), so the line is matched as a whole rather than split.
LOG_ROW = re.compile(
    r"^\s*S(\d+)\s*/\s*T(\d+)\(Ms=\s*([+-]?\d+)\)\s+"
    r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
)
LOG_D4 = re.compile(r"^\s*PyOQP dftd correction (-?\d+\.\d+)\s*$", re.MULTILINE)
#: The SCF energy before any dispersion correction, as the manifold section
#: prints it.  SOC always runs on a triplet ROHF reference, so this is the one
#: starting point both manifolds share.  The line recurs once uncorrected and
#: once corrected, so the FIRST match is the uncorrected one.
LOG_REFERENCE = re.compile(
    r"^\s*PyOQP triplet ROHF reference \(internal\)\s+(-?\d+\.\d+)\s*$",
    re.MULTILINE)


def _oqp_available():
    try:
        from oqp import lib  # noqa: F401
        return hasattr(lib, "soc_mrsf")
    except Exception:
        return False


def _parse_log_rows(log_text):
    rows = []
    for line in log_text.splitlines():
        match = LOG_ROW.match(line)
        if match:
            ist, jst, ms, re1, im1, re2, im2, absolute = match.groups()
            rows.append({
                's': int(ist), 't': int(jst), 'ms': int(ms),
                're': float(re1) + float(re2), 'im': float(im1) + float(im2),
                'absolute': float(absolute),
            })
    return rows


def _run_deck_text(text, workdir, log_name):
    """Run one deck inside ``workdir`` and return the path of its log."""
    from oqp.pyoqp import Runner

    inp = Path(workdir) / f'{log_name}.inp'
    inp.write_text(text, encoding='utf-8')
    log = Path(workdir) / f'{log_name}.log'
    Runner(input_file=str(inp), log=str(log), silent=1, usempi=False).run()
    return log


def _run_deck(deck, workdir, log_name='run'):
    return _run_deck_text(Path(deck).read_text(encoding='utf-8'), workdir, log_name)


@unittest.skipUnless(_oqp_available(), "OpenQP shared library is not built")
class SOCExportTests(unittest.TestCase):
    """Everything the SOC step of runtype=prop writes."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory(prefix='soc_export_')
        cls.workdir = Path(cls._tmp.name)
        cls.log = _run_deck(SOC_DECK, cls.workdir, log_name='soc_prop')
        cls.log_text = cls.log.read_text(encoding='utf-8')

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def _table(self):
        path = self.workdir / 'socme'
        self.assertTrue(path.is_file(), 'export=true must write socme')
        return np.loadtxt(path, usecols=(0, 1, 2, 3, 4, 5), ndmin=2)

    def test_socme_reproduces_the_native_singlet_triplet_table(self):
        table = self._table()
        rows = _parse_log_rows(self.log_text)
        self.assertGreater(len(rows), 0, 'no S-T table found in the log')
        self.assertEqual(len(table), len(rows),
                         f'socme has {len(table)} rows, the log {len(rows)}')

        # The log prints each of the four components at 4 decimals, so their sum
        # carries up to 1e-4 of formatting error; Absolute is printed at 12.
        for index, (line, expected) in enumerate(zip(table, rows)):
            ist, jst, ms, value_re, value_im, absolute = line
            self.assertEqual((int(ist), int(jst), int(ms)),
                             (expected['s'], expected['t'], expected['ms']),
                             f'row {index}: labels differ from the log')
            self.assertLess(abs(value_re - expected['re']), 1.1e-4,
                            f'row {index}: Re {value_re} vs log {expected["re"]}')
            self.assertLess(abs(value_im - expected['im']), 1.1e-4,
                            f'row {index}: Im {value_im} vs log {expected["im"]}')
            self.assertLess(abs(absolute - expected['absolute']), 1e-6,
                            f'row {index}: Absolute {absolute} vs log '
                            f'{expected["absolute"]}')
            self.assertLess(abs(absolute - np.hypot(value_re, value_im)), 1e-9,
                            f'row {index}: Absolute is not |Re + i Im|')

    def test_socme_carries_imaginary_couplings_so_the_sign_check_is_live(self):
        # Reading the two-dimensional tag without transposing it back keeps the
        # Fortran matrix but flips every imaginary entry.  That is only caught
        # above if the imaginary part is actually resolved by the tolerance.
        table = self._table()
        self.assertGreater(np.abs(table[:, 4]).max(), 0.1,
                           'socme Im column is numerically zero: the comparison '
                           'against the log cannot detect a sign flip')

    def test_ms_triple_is_complete_for_every_state_pair(self):
        table = self._table()
        labels = {(int(row[0]), int(row[1]), int(row[2])) for row in table}
        nstate_s = max(int(row[0]) for row in table) + 1
        nstate_t = max(int(row[1]) for row in table) + 1
        self.assertEqual(
            labels,
            {(ist, jst, ms) for ist in range(nstate_s)
             for jst in range(nstate_t) for ms in (-1, 0, 1)},
            'socme must hold one row per (singlet, triplet, Ms) triple')

    def test_each_spin_manifold_is_published_once(self):
        singlet = np.loadtxt(self.workdir / 'energies_s')
        triplet = np.loadtxt(self.workdir / 'energies_t')
        self.assertEqual(singlet.shape, triplet.shape)
        self.assertFalse(np.allclose(singlet[1:], triplet[1:]),
                         'the two ladders must differ, or the labels were mixed')

        corrections = self.d4_corrections()
        self.assertGreater(abs(corrections[0]), 1e-7,
                           'the deck must carry a non-zero dispersion '
                           'correction for this to check anything')
        match = LOG_REFERENCE.search(self.log_text)
        self.assertIsNotNone(match, 'no uncorrected reference energy in the log')
        scf = float(match.group(1))
        # LastStep adds the dispersion in place, so a manifold that passed
        # through it twice would land a whole correction below this line.
        for manifold, energies in (('singlet', singlet), ('triplet', triplet)):
            self.assertAlmostEqual(
                energies[0], scf + corrections[0], delta=1e-8,
                msg=f'{manifold} ground state is not SCF plus one dispersion '
                    'correction')

        self.assertFalse(np.allclose(np.loadtxt(self.workdir / 'grad_s_2'),
                                     np.loadtxt(self.workdir / 'grad_t_2')),
                         'the two manifolds published the same gradient')

    def d4_corrections(self):
        """The dispersion corrections the log reports, in Hartree."""
        values = [float(match) for match in LOG_D4.findall(self.log_text)]
        self.assertGreaterEqual(len(values), 2,
                                'expected one dispersion line per manifold')
        return values


@unittest.skipUnless(_oqp_available(), "OpenQP shared library is not built")
class SOCStandaloneExportTests(unittest.TestCase):
    """runtype=soc gains the same export, and exercises two-digit labels."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory(prefix='soc_standalone_')
        cls.workdir = Path(cls._tmp.name)
        deck = SOC_DECK_STANDALONE.read_text(encoding='utf-8').rstrip()
        cls.log = _run_deck_text(deck + '\n\n[properties]\nexport=true\n',
                                 cls.workdir, 'standalone')
        cls.log_text = cls.log.read_text(encoding='utf-8')

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_standalone_soc_takes_no_gradients(self):
        """runtype=soc solves both manifolds and couples them; nothing gradients."""
        published = sorted(p.name for p in self.workdir.iterdir()
                           if p.name.startswith('grad'))
        self.assertEqual(published, [],
                         f'runtype=soc wrote gradient files: {published}')

    def test_standalone_soc_exports_the_couplings_and_both_ladders(self):
        here = self.workdir
        for name in ('socme', 'energies_s', 'energies_t'):
            self.assertTrue((here / name).is_file(), f'{name} missing')

        table = np.loadtxt(here / 'socme', ndmin=2)
        rows = _parse_log_rows(self.log_text)
        # ns = nt = 12 here, so the singlet index reaches two digits: the label
        # encoding of the file has to agree with the log beyond single digits.
        self.assertEqual(len(table), 12 * 12 * 3)
        self.assertEqual(len(rows), len(table))
        self.assertIn(10, {int(row[0]) for row in table},
                      'no two-digit singlet index was exercised')
        for line, expected in zip(table, rows):
            self.assertEqual((int(line[0]), int(line[1]), int(line[2])),
                             (expected['s'], expected['t'], expected['ms']))
            # Each component is printed at four decimals, so two rounded parts
            # sum with up to 1e-4 of formatting error; Absolute is at 12.
            self.assertLess(abs(line[3] - expected['re']), 1.1e-4)
            self.assertLess(abs(line[4] - expected['im']), 1.1e-4)
            self.assertLess(abs(line[5] - expected['absolute']), 1e-6)


GEOMETRY_PREVIOUS = """system=
 8   0.000000000   0.000000000   0.158938446
 1  -0.633194329   0.633194329  -0.514469223
 1   0.633194329  -0.633194329  -0.514469223
"""
SYSTEM_BLOCK = re.compile(r'^system=\n(?: .*\n){3}', re.MULTILINE)


class SOCManifoldCouplingTests(unittest.TestCase):
    """Each manifold couples against its OWN previous-step amplitudes.

    Two consecutive spin-orbit steps wired the way an external driver wires an
    MD step (mol.back_door carries the earlier geometry and its data), plus a run
    that only ever solved the triplet manifold as the oracle for those couplings.

    This deck keeps its two response roots far apart, so the true coupling is
    ~1e-12 and the fixture is not judged by magnitude: pairing the triplet with
    the singlet history instead was measured to move that entry to 3.2e-2, four
    orders above the oracle, so agreeing with the oracle to 1e-6 only holds if
    each manifold used its own history.
    """

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory(prefix='soc_nacme_')
        cls.workdir = Path(cls._tmp.name)
        base = SOC_DECK.read_text(encoding='utf-8')
        # Unequal ladder lengths, so the two published matrices cannot be the
        # same array pretending to be two manifolds.
        soc_deck = (base.replace('nstate=2\n', 'nstate=2\nnstate_s=2\nnstate_t=3\n')
                    .replace('grad_t=2', 'grad_t=3'))
        triplet_only = (base.replace('soc=true\n', '')
                        .replace('grad_s=2\ngrad_t=2', 'grad=3')
                        .replace('[tdhf]\ntype=mrsf\nnstate=2\n',
                                 '[tdhf]\ntype=mrsf\nnstate=3\nmultiplicity=3\n'))

        # The earlier step has no history yet, so it publishes the checkpoint and
        # computes no couplings.
        first = SYSTEM_BLOCK.sub(GEOMETRY_PREVIOUS, soc_deck, count=1)
        first = first.replace('nac=nacme\n', '')
        folder = cls.workdir / 'step1'
        folder.mkdir()
        cls.previous = cls._run(first, folder, 'step1')
        cls.soc = cls._run(soc_deck, cls.workdir, 'step2', cls.previous)

        # The oracle is the triplet manifold from a step that only ever solved
        # triplets, so its checkpoint's generic amplitudes are that manifold's.
        # A singlet oracle would need a singlet-only previous step as well; the
        # singlet block is the pre-existing path, covered by the NACME tests.
        reference = cls.workdir / 'ref_triplet'
        reference.mkdir()
        cls._run(triplet_only, reference, 'ref', cls.previous)
        cls.oracle = np.loadtxt(reference / 'nacme')

    @classmethod
    def _run(cls, text, folder, name, history=None):
        from oqp.pyoqp import Runner

        deck = folder / f'{name}.inp'
        deck.write_text(text, encoding='utf-8')
        runner = Runner(input_file=str(deck), log=str(folder / f'{name}.log'),
                        silent=1, usempi=False)
        if history is not None:
            runner.mol.back_door = history
        runner.run()
        return (np.array(runner.mol.get_system(), copy=True),
                runner.mol.get_data())

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_both_manifolds_publish_their_own_coupling_matrices(self):
        for name in ('nacme_s', 'nacme_t', 'dcme_s', 'dcme_t'):
            self.assertTrue((self.workdir / name).is_file(), f'{name} missing')

        produced = {'s': np.loadtxt(self.workdir / 'nacme_s'),
                    't': np.loadtxt(self.workdir / 'nacme_t')}
        # The two ladders have different lengths here, so one array cannot stand
        # in for both.
        self.assertEqual(produced['s'].shape, (2, 2))
        self.assertEqual(produced['t'].shape, (3, 3))
        difference = np.abs(produced['t'] - self.oracle).max()
        self.assertLess(difference, 1e-6,
                        'nacme_t differs from the manifold-only triplet run by '
                        f'{difference:.3e}')

    def test_the_history_really_is_a_different_step(self):
        """Guards against the agreement coming from pairing a step with itself."""
        current = np.asarray(self.soc[1]['OQP::td_bvec_mo_t'])
        previous = np.asarray(self.previous[1]['OQP::td_bvec_mo_t'])
        self.assertEqual(current.shape, previous.shape)
        self.assertGreater(np.abs(current - previous).max(), 1e-3,
                           'the two steps produced the same triplet amplitudes')


class TwoLadderHistoryWithoutSOCTests(unittest.TestCase):
    """A two-ladder step without spin-orbit coupling still records its history.

    The unlabelled amplitude tags end up holding the ladder solved last (the
    triplet), so the next step's singlet NACME must read the labelled singlet
    record.  Without one, it overlapped singlet amplitudes with the previous
    step's triplets, measured 1.3e-1 off the singlet-only oracle here with
    nothing reported.  Agreeing with that oracle to 1e-6 therefore only holds
    if the singlets were paired with their own history.
    """

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory(prefix='history_nosoc_')
        cls.workdir = Path(cls._tmp.name)
        deck = (SOC_DECK.read_text(encoding='utf-8')
                .replace('soc=true\n', '')
                .replace('basis=6-31G*', 'basis=STO-3G')
                .replace('d4=true\n', ''))
        single = deck.replace('grad_t=2\n', '').replace('grad_s=2', 'grad=2')
        run = SOCManifoldCouplingTests._run
        for name, text in (('both', deck), ('singlet', single)):
            folder = cls.workdir / name
            (folder / 'step1').mkdir(parents=True)
            first = SYSTEM_BLOCK.sub(GEOMETRY_PREVIOUS, text, count=1)
            previous = run(first.replace('nac=nacme\n', ''), folder / 'step1', 'step1')
            run(text, folder, 'step2', previous)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_singlet_couplings_use_the_singlet_history(self):
        produced = np.loadtxt(self.workdir / 'both' / 'nacme_s', ndmin=2)
        oracle = np.loadtxt(self.workdir / 'singlet' / 'nacme', ndmin=2)
        self.assertEqual(produced.shape, oracle.shape)
        difference = np.abs(produced - oracle).max()
        self.assertLess(difference, 1e-6,
                        'nacme_s differs from the singlet-only run by '
                        f'{difference:.3e}')


PER_LADDER_ROOTS_DECK = """# per-ladder root counts only: nstate_s / nstate_t, no nstate
[input]
system=
 8   0.000000000   0.000000000  -0.041061554
 1  -0.533194329   0.533194329  -0.614469223
 1   0.533194329  -0.533194329  -0.614469223
runtype=prop
method=tdhf
functional=bhhlyp
basis=6-31G*
charge=0
ispher=false

[scf]
type=rohf
multiplicity=3
converger_type=soscf
maxit=200
conv=1e-8

[tdhf]
type=mrsf
nstate_s=3
nstate_t=4
maxit=100
maxit_zv=100
zvconv=1.0e-10

[guess]
type=huckel
save_mol=True

[properties]
export=true
soc=true
back_door=true
grad_s=2
grad_t=2
"""


class SOCRootCountTests(unittest.TestCase):
    """nstate_s / nstate_t size the ladder being answered, both rounds.

    Before the fix only the second round was sized: the round that answered the
    configured multiplicity ran on [tdhf] nstate, whose default is 1, so a deck
    naming per-ladder root counts converged one root and then refused the
    gradient it had been asked for.
    """

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory(prefix='roots_')
        cls.workdir = Path(cls._tmp.name)
        _run_deck_text(PER_LADDER_ROOTS_DECK, cls.workdir, 'roots')

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_each_ladder_converged_the_number_of_roots_it_named(self):
        energies_s = np.loadtxt(self.workdir / 'energies_s')
        energies_t = np.loadtxt(self.workdir / 'energies_t')
        # one reference row plus the roots of that ladder
        self.assertEqual(energies_s.shape, (1 + 3,))
        self.assertEqual(energies_t.shape, (1 + 4,))
        # the couplings cover every (singlet, triplet, Ms) triple
        table = np.loadtxt(self.workdir / 'socme', ndmin=2)
        self.assertEqual(table.shape[0], 3 * 4 * 3)
        self.assertTrue((table[:, 0] <= 2).all() and (table[:, 1] <= 3).all())

    def test_the_configured_ladder_gradient_was_taken(self):
        singlet = np.loadtxt(self.workdir / 'grad_s_2')
        triplet = np.loadtxt(self.workdir / 'grad_t_2')
        self.assertEqual(singlet.shape, triplet.shape)
        self.assertTrue(np.isfinite(singlet).all() and np.isfinite(triplet).all())
        self.assertFalse(np.allclose(singlet, triplet),
                         'both ladders published one gradient')

    def test_a_shorter_ladder_is_refused_before_the_gradient(self):
        # The witness that the run above is not a fixture artifact: ask for root 2
        # of a ladder that converges 1 root and the request is refused.
        from oqp.utils.input_checker import check_input_values

        config = {
            'input': {'runtype': 'prop', 'method': 'tdhf', 'functional': 'bhhlyp',
                      'system': '\nO 0.0 0.0 0.0\nH 0.0 0.0 0.95\nH 0.9 0.0 -0.3'},
            'scf': {'type': 'rohf', 'multiplicity': 3},
            'tdhf': {'type': 'mrsf', 'nstate_s': 1, 'nstate_t': 4},
            'properties': {'soc': True, 'grad_s': [2]},
        }
        report = check_input_values(config, raise_error=False, emit=False)
        self.assertIn('tdhf.nstate_s', {item.path for item in report.errors},
                      report.to_text())


class SOCBranchTests(unittest.TestCase):
    """The paths a single deck cannot cover."""

    def test_gradient_run_publishes_both_ladders(self):
        """runtype=grad honours grad_s and grad_t and labels what it writes."""
        deck = (SOC_DECK.read_text(encoding='utf-8')
                .replace('runtype=prop', 'runtype=grad')
                .replace('soc=true\n', '')
                .replace('nac=nacme\n', '')
                .replace('back_door=true\n', ''))
        with tempfile.TemporaryDirectory(prefix='grad_both_') as workdir:
            _run_deck_text(deck, workdir, 'grad')
            here = Path(workdir)
            for name in ('energies_s', 'energies_t', 'grad_s_2', 'grad_t_2'):
                self.assertTrue((here / name).is_file(), f'{name} missing')
            self.assertFalse((here / 'socme').is_file(),
                             'runtype=grad computes no spin-orbit coupling')
            singlet = np.loadtxt(here / 'grad_s_2')
            triplet = np.loadtxt(here / 'grad_t_2')
            self.assertFalse(np.allclose(singlet, triplet),
                             'both labels published one manifold gradient')

    def test_a_second_ladder_is_published_without_spin_orbit(self):
        """grad_t alone solves and publishes that ladder, and asks for no socme."""
        deck = (PLAIN_PROP_DECK.read_text(encoding='utf-8')
                .replace('grad=2', 'grad=2\ngrad_t=1'))
        with tempfile.TemporaryDirectory(prefix='soc_ladder_') as workdir:
            _run_deck_text(deck, workdir, 'ladder')
            here = Path(workdir)
            for name in ('energies_s', 'energies_t', 'grad_s_2', 'grad_t_1',
                         'nacme_s', 'dcme_s'):
                self.assertTrue((here / name).is_file(), f'{name} missing')
            self.assertFalse((here / 'socme').is_file(),
                             'no spin-orbit coupling was requested')
            singlet = np.loadtxt(here / 'energies_s')
            triplet = np.loadtxt(here / 'energies_t')
            self.assertEqual(singlet.shape, triplet.shape)
            self.assertFalse(np.allclose(singlet[1:], triplet[1:]),
                             'the second ladder duplicated the first one')

    def test_run_without_soc_keeps_the_unlabelled_names(self):
        # The unlabelled files are the engine's own channel between a numerical
        # Hessian or NAC parent and its displacement workers, so a properties
        # run that asks for no spin-orbit coupling must not rename them.
        labelled = ('energies_s', 'energies_t', 'grad_s_2', 'grad_t_2',
                    'nacme_s', 'nacme_t', 'dcme_s', 'dcme_t', 'socme')
        with tempfile.TemporaryDirectory(prefix='soc_plain_') as workdir:
            _run_deck(PLAIN_PROP_DECK, workdir, log_name='plain')
            here = Path(workdir)
            for name in ('energies', 'grad_2', 'nacme', 'dcme'):
                self.assertTrue((here / name).is_file(), f'{name} disappeared')
            present = {path.name for path in here.iterdir()}
            self.assertEqual(present & set(labelled), set(),
                             f'a run without properties.soc wrote {labelled}')

    def test_a_triplet_only_request_does_not_differentiate_the_reference(self):
        # soc=true with grad_t alone: the singlet ladder has no request.  The
        # [properties] grad default (root 0, the internal ROHF reference) used
        # to fill it, so the run took a Z-vector of the reference, which does
        # not converge, and stopped before publishing anything.
        deck = (SOC_DECK.read_text(encoding='utf-8')
                .replace('runtype=prop', 'runtype=data')
                .replace('basis=6-31G*', 'basis=STO-3G')
                .replace('nac=nacme\n', '')
                .replace('grad_s=2\n', '')
                .replace('grad_t=2', 'grad_t=1')
                .replace('d4=true\n', ''))
        with tempfile.TemporaryDirectory(prefix='soc_triplet_only_') as workdir:
            log = _run_deck_text(deck, workdir, 'triplet_only')
            here = Path(workdir)
            for name in ('socme', 'energies_s', 'energies_t', 'grad_t_1'):
                self.assertTrue((here / name).is_file(), f'{name} missing')
            self.assertEqual(sorted(path.name for path in here.glob('grad_s_*')), [],
                             'no singlet gradient was requested')
            self.assertNotIn('Gradient of triplet ROHF reference',
                             log.read_text(encoding='utf-8'))

    def test_gradient_run_sizes_the_configured_ladder_from_its_own_count(self):
        # runtype=grad answered its configured ladder on [tdhf] nstate even when
        # nstate_s named another count, while the checker (like prop and data)
        # judged grad_s against nstate_s: a deck the checker passed then refused
        # the gradient after the solve.
        deck = (SOC_DECK.read_text(encoding='utf-8')
                .replace('runtype=prop', 'runtype=grad')
                .replace('basis=6-31G*', 'basis=STO-3G')
                .replace('soc=true\n', '')
                .replace('nac=nacme\n', '')
                .replace('back_door=true\n', '')
                .replace('nstate=2', 'nstate=1\nnstate_s=2')
                .replace('grad_t=2', 'grad_t=1')
                .replace('d4=true\n', ''))
        with tempfile.TemporaryDirectory(prefix='grad_counts_') as workdir:
            _run_deck_text(deck, workdir, 'counts')
            here = Path(workdir)
            # one reference row plus the roots of that ladder
            self.assertEqual(np.loadtxt(here / 'energies_s').shape, (1 + 2,))
            self.assertEqual(np.loadtxt(here / 'energies_t').shape, (1 + 1,))
            for name in ('grad_s_2', 'grad_t_1'):
                self.assertTrue((here / name).is_file(), f'{name} missing')

    def test_data_runtype_publishes_the_second_manifold_and_socme(self):
        # The same request on runtype=data, one basis smaller and without a NAC
        # section, because runtype=data takes no nacme branch.
        deck = (SOC_DECK.read_text(encoding='utf-8')
                .replace('runtype=prop', 'runtype=data')
                .replace('basis=6-31G*', 'basis=STO-3G')
                .replace('nac=nacme\n', '')
                .replace('grad_s=2', 'grad_s=1')
                .replace('grad_t=2', 'grad_t=1')
                .replace('d4=true\n', ''))
        with tempfile.TemporaryDirectory(prefix='soc_data_') as workdir:
            deck_path = Path(workdir) / 'data.inp'
            deck_path.write_text(deck, encoding='utf-8')
            _run_deck(deck_path, workdir, log_name='data')
            here = Path(workdir)
            for name in ('socme', 'energies_s', 'energies_t',
                         'grad_s_1', 'grad_t_1'):
                self.assertTrue((here / name).is_file(), f'{name} missing')
            rows = _parse_log_rows((here / 'data.log').read_text(encoding='utf-8'))
            table = np.loadtxt(here / 'socme', usecols=(3, 4, 5), ndmin=2)
            self.assertEqual(len(table), len(rows))
            self.assertLess(abs(table[:, 2] - [row['absolute'] for row in rows]).max(),
                            1e-6)


class SOCFileFormatTests(unittest.TestCase):
    """Naming of the new export, without a calculation behind it."""

    def test_socme_name_follows_the_export_title(self):
        from oqp.utils.file_utils import dump_data

        class Mol:
            usempi = False

        table = np.array([[0, 0, -1, 1.0, 2.0, np.hypot(1.0, 2.0)]])
        with tempfile.TemporaryDirectory(prefix='socme_io_') as workdir:
            dump_data(Mol(), (table, 'sys1'), title='SOCME', fpath=workdir)
            dump_data(Mol(), (table, ''), title='SOCME', fpath=workdir)
            self.assertTrue((Path(workdir) / 'sys1.socme').is_file())
            self.assertTrue((Path(workdir) / 'socme').is_file())
            back = np.loadtxt(Path(workdir) / 'socme', ndmin=2)
            self.assertEqual(back.shape, (1, 6))

    def test_spin_labels_never_collide_with_the_unlabelled_names(self):
        from oqp.utils.file_utils import _export_stem

        self.assertEqual(_export_stem('ENERGY'), 'energies')
        self.assertEqual(_export_stem('GRADIENT'), 'grad')
        self.assertEqual(_export_stem('ENERGY', 's'), 'energies_s')
        self.assertEqual(_export_stem('GRADIENT', 't'), 'grad_t')
        self.assertEqual(_export_stem('NACV', 's'), 'nac_s')
        self.assertEqual(_export_stem('SOCME'), 'socme')
        with self.assertRaises(ValueError):
            _export_stem('ENERGY', 'q')


class SOCRankGuardTests(unittest.TestCase):
    """Only the world root holds the SOC tags, and it must be the only reader."""

    class _NoTags:
        """Stand-in for a rank that never received OQP::soc_*."""

        def __init__(self):
            self.reads = []

        def __getitem__(self, key):
            self.reads.append(key)
            raise AttributeError(f'Key `{key}` not found in QOP data')

        def set_tdhf_multiplicity(self, *_):
            pass

        def set_tdhf_nstate(self, *_):
            pass

    def _soc(self, world_rank):
        from unittest import mock

        import oqp.library.single_point as single_point

        config = {
            'input': {'method': 'tdhf'},
            'tdhf': {'nstate': 2, 'nstate_s': 0, 'nstate_t': 0,
                     'multiplicity': 1, 'type': 'mrsf'},
            'properties': {'export': True, 'title': ''},
            'guess': {'save_mol': False},
            'tests': {'exception': False},
        }
        mol = type('Mol', (), {'config': config, 'data': self._NoTags(),
                               'usempi': True})()
        with mock.patch.object(single_point, 'MPIManager',
                               return_value=type('Manager', (),
                                                 {'world_rank': world_rank})()):
            soc = single_point.SOC(mol)
            soc.export_socme()
        return mol.data.reads

    def test_a_non_root_rank_never_reads_the_soc_tags(self):
        # The native kernel allocates OQP::soc_* only on its own rank, and tag
        # arrays are per rank, so a read from another rank raises.  @mpi_dump on
        # dump_data guards the write, not the argument building in front of it.
        self.assertEqual(self._soc(1), [])

    def test_the_root_rank_does_read_them(self):
        # Without this half the guard above would pass for a method that simply
        # never exports anything.
        with self.assertRaises(AttributeError):
            self._soc(0)


class SOCConstantTests(unittest.TestCase):
    """One conversion factor, three places that compute it."""

    def test_dfac_matches_the_tight_binding_and_namd_sites(self):
        from oqp.utils.constants import (FINE_STRUCTURE, HA_TO_WAVENUM,
                                         SOC_AU_TO_WAVENUM)
        from oqp.library import namd, openqp_dftb

        expected = 0.5 * FINE_STRUCTURE ** 2 * HA_TO_WAVENUM
        self.assertEqual(SOC_AU_TO_WAVENUM, expected)
        self.assertEqual(HA_TO_WAVENUM, openqp_dftb.HA_TO_WAVENUMBER)
        self.assertEqual(HA_TO_WAVENUM, namd.HA_TO_WAVENUM)
        self.assertEqual(FINE_STRUCTURE, openqp_dftb.FINE_STRUCTURE)


class SOCGateTests(unittest.TestCase):
    """What properties.soc refuses, before any calculation is paid for."""

    @staticmethod
    def _config(runtype='prop'):
        return {
            'input': {
                'runtype': runtype,
                'method': 'tdhf',
                'functional': 'bhhlyp',
                'basis': '6-31g*',
                'system': '\nO 0.0 0.0 0.0\nH 0.0 0.0 0.95\nH 0.9 0.0 -0.3',
            },
            'guess': {},
            'scf': {'type': 'rohf', 'multiplicity': 3},
            'tdhf': {'type': 'mrsf', 'nstate': 2},
            'properties': {'soc': True, 'grad_s': [1], 'grad_t': [1]},
            'nac': {'states': [[1, 2]]},
            'optimize': {},
        }

    @staticmethod
    def _report(config):
        from oqp.utils.input_checker import check_input_values
        return check_input_values(config, raise_error=False, emit=False)

    def _paths(self, config):
        return {item.path for item in self._report(config).errors}

    def test_a_complete_request_passes(self):
        self.assertEqual(self._paths(self._config()), set(),
                         self._report(self._config()).to_text())
        self.assertEqual(self._paths(self._config('data')), set())

    def test_runtype_that_never_reaches_the_soc_step_is_refused(self):
        self.assertIn('properties.soc', self._paths(self._config('energy')))

    def test_the_reference_requirements_are_named(self):
        config = self._config()
        config['scf']['type'] = 'rhf'
        config['scf']['multiplicity'] = 1
        config['tdhf']['type'] = 'sf'
        errors = {item.path: item.message
                  for item in self._report(config).errors}
        self.assertIn('scf.type', errors)
        self.assertIn('scf.multiplicity', errors)
        self.assertIn('tdhf.type', errors)
        self.assertIn('properties.soc', ' '.join(errors.values()),
                      'the message must say which switch asked for SOC')

    def test_core_shift_is_refused(self):
        config = self._config()
        config['tdhf']['ixcore'] = '0'
        self.assertIn('tdhf.ixcore', self._paths(config))

    def test_two_electron_soc_with_spherical_ao_is_refused_early(self):
        # cc-pVDZ publishes spherical d shells; 6-31g* publishes Cartesian ones.
        for basis, flagged in (('cc-pvdz', True), ('6-31g*', False)):
            config = self._config()
            config['input']['basis'] = basis
            if flagged:
                self.assertIn('input.soc_2e', self._paths(config), basis)
            else:
                self.assertEqual(self._paths(config), set(), basis)
            # Dropping the two-electron term is the documented way out.
            config['input']['soc_2e'] = 0
            self.assertNotIn('input.soc_2e', self._paths(config), basis)

    def test_gradient_beyond_the_manifold_state_count_is_refused(self):
        config = self._config()
        config['tdhf']['nstate_t'] = 1
        config['properties']['grad_t'] = [2]
        self.assertIn('tdhf.nstate_t', self._paths(config))
        # The singlet list is unaffected by the triplet's shortfall.
        config['tdhf']['nstate_s'] = 2
        config['properties']['grad_s'] = [2]
        self.assertNotIn('tdhf.nstate_s', self._paths(config))

    def test_a_spin_orbit_run_does_not_read_the_shared_grad_selector(self):
        config = self._config()
        config['properties']['grad'] = [2]
        self.assertIn('properties.grad', self._paths(config))
        # The schema default (root 0) is not an explicit request.
        config['properties']['grad'] = [0]
        self.assertEqual(self._paths(config), set(),
                         self._report(config).to_text())

    def test_nac_vectors_are_not_offered_per_manifold(self):
        # The analytic engine has no triplet fold (nac_analytic raises
        # NotImplementedError), so a spin-orbit run must not promise vectors.
        for route in ('numerical', 'analytical'):
            config = self._config()
            config['properties']['nac'] = 'nac'
            config['nac'] = {'states': [[1, 2]], 'type': route}
            self.assertIn('properties.nac', self._paths(config), route)
        config = self._config()
        config['properties']['nac'] = 'nacme'
        self.assertNotIn('properties.nac', self._paths(config))

    def test_runtype_soc_refuses_gradient_selectors(self):
        for key in ('grad_s', 'grad_t'):
            config = self._config('soc')
            config['properties'] = {key: [2]}
            self.assertIn(f'properties.{key}', self._paths(config), key)
        config = self._config('soc')
        config['properties'] = {'grad_t': [2]}
        config['input']['runtype'] = 'data'
        self.assertNotIn('properties.grad_t', self._paths(config))

    def test_a_second_ladder_can_be_asked_for_without_spin_orbit(self):
        # grad_t names a ladder outright, so it is not a spin-orbit-only knob.
        config = self._config()
        config['properties'] = {'grad': [0], 'grad_t': [1]}
        self.assertEqual(self._paths(config), set(),
                         self._report(config).to_text())

    def test_a_spin_orbit_run_never_falls_back_to_the_shared_selector(self):
        # The schema default grad=[0] names the internal reference, not a state:
        # a spin-orbit run with grad_t alone asks for no singlet gradient.
        from oqp.utils.spin_states import grad_requests

        config = self._config()
        config['properties'] = {'soc': True, 'grad': [0], 'grad_t': [1]}
        self.assertEqual(grad_requests(config), {'s': [], 't': [1]})
        self.assertEqual(self._paths(config), set(),
                         self._report(config).to_text())
        # Without spin-orbit coupling the shared selector keeps describing the
        # configured ladder, as every existing deck relies on.
        config['properties'] = {'grad': [2], 'grad_t': [1]}
        self.assertEqual(grad_requests(config), {'s': [2], 't': [1]})

    def test_tight_binding_spin_orbit_is_only_offered_standalone(self):
        # The properties/data flow switches ladders through [tdhf]
        # multiplicity, which the tight-binding adapters do not read.
        for method in ('dftb', 'xtb'):
            for runtype in ('prop', 'data'):
                config = self._config(runtype)
                config['input']['method'] = method
                config['properties'] = {'soc': True}
                self.assertIn('properties.soc', self._paths(config),
                              f'{method} {runtype}')
            config = self._config('soc')
            config['input']['method'] = method
            config['properties'] = {}
            self.assertNotIn('properties.soc', self._paths(config), method)

    def test_per_manifold_selectors_need_the_mrsf_ladders(self):
        # Without spin-orbit coupling as well: the second ladder is reached by
        # switching [tdhf] multiplicity, which only the MRSF response answers.
        for method, td_type in (('dftb', 'mrsf'), ('xtb', 'mrsf'), ('tdhf', 'rpa'),
                                ('hf', 'mrsf')):
            config = self._config('prop')
            config['input']['method'] = method
            config['tdhf']['type'] = td_type
            config['properties'] = {'grad_t': [1]}
            self.assertIn('properties.grad_t', self._paths(config),
                          f'{method}/{td_type}')
        # A runtype that never reads the selectors must not drop them quietly.
        config = self._config('energy')
        config['properties'] = {'grad_t': [1]}
        self.assertIn('properties.grad_t', self._paths(config))
        # The root count is checked without spin-orbit coupling too.
        config = self._config('grad')
        config['tdhf']['nstate_s'] = 1
        config['properties'] = {'grad_s': [2], 'grad_t': [1]}
        self.assertIn('tdhf.nstate_s', self._paths(config))
        # Roots count from 1; 0 and below would reach the native Z-vector.
        for states in ([0], [-1, 2]):
            config = self._config()
            config['properties']['grad_t'] = states
            self.assertIn('properties.grad_t', self._paths(config), states)

    def test_data_runtype_refuses_nacme(self):
        # runtype=data has no overlap with the previous step, so a nacme
        # request there would finish without its coupling matrices.
        config = self._config('data')
        config['properties']['nac'] = 'nacme'
        self.assertIn('properties.nac', self._paths(config))
        config = self._config('prop')
        config['properties']['nac'] = 'nacme'
        self.assertNotIn('properties.nac', self._paths(config))


class SOCExampleReferenceTests(unittest.TestCase):
    """The reference keys a lean dump is allowed to drop must stay."""

    def test_standalone_soc_references_still_publish_soc(self):
        import json

        # regression.py cannot mark 'soc' required for every runtype that now
        # carries it, so the standalone references are pinned here instead.
        for name in ('H2O_BHHLYP_SOC', 'CH3Br-BHHLYP-SOC',
                     'H2O_BHHLYP_MRSF_PROP_SOC'):
            reference = json.loads(
                (EXAMPLES / 'SOC' / f'{name}.json').read_text(encoding='utf-8'))
            self.assertIn('soc', reference, name)
            self.assertIn('td_singlet_energies', reference, name)
            self.assertIn('td_triplet_energies', reference, name)

    def test_the_example_turns_the_flag_on(self):
        # openqp --check_feature_coverage only reads .inp files.
        text = SOC_DECK.read_text(encoding='utf-8')
        self.assertRegex(text, r'(?m)^\s*soc\s*=\s*true\s*$')


class SOCPythonAPITests(unittest.TestCase):
    """Rule 3: the keyword is reachable from the compact API."""

    def test_properties_soc_keyword_is_exposed(self):
        from oqp.openqp import OpenQP

        job = OpenQP(project='soc_api')
        job.properties(soc=True)
        self.assertTrue(job.to_input_dict()['properties']['soc'])

        job = OpenQP(project='soc_api_settings')
        job.settings.properties(soc=True, export=True)
        config = job.to_input_dict()
        self.assertTrue(config['properties']['soc'])
        self.assertTrue(config['properties']['export'])

    def test_concise_input_spelling_lowers_to_the_keyword(self):
        from oqp.utils import oqp_input

        deck = EXAMPLES / 'SOC' / 'H2O_BHHLYP_MRSF_PROP_SOC.oqp'
        spec = oqp_input.parse_canonical_oqp(deck.read_text(encoding='utf-8'))
        config = oqp_input.lower_to_legacy(spec, source_dir=deck.parent)
        self.assertEqual(config['input']['runtype'], 'prop')
        self.assertEqual(config['properties']['soc'].lower(), 'true')
        self.assertEqual(config['properties']['export'].lower(), 'true')
        # prop(S1) names the singlet ladder, so it lowers to grad_s and never to
        # the ambiguous shared selector.
        self.assertEqual(config['properties']['grad_s'], '2')
        self.assertNotIn('grad', config['properties'])


if __name__ == '__main__':
    unittest.main()
