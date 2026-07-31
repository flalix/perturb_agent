"""Testes dos guards. Se algum destes quebrar, o sistema voltou a ser um
gerador de narrativas plausíveis. Rode em CI."""

import math
import unittest

from sciagent_core import (
    Actor, Assay, Claim, Direction, EntityType, EvidenceItem, EvidenceKind,
    GuardError, Hypothesis, Ledger, Prediction, RefutationAttempt, Scope,
    SeedResolver, State, Triple, attenuate_lr, expected_information_gain,
)
from sciagent_core.evidence import DependenceModel, aggregate_log_lr, information_ceiling

R = SeedResolver()
GEN = Actor("g", "generator", "A", "1")
CRIT1 = Actor("c1", "critic", "B", "1")
CRIT2 = Actor("c2", "critic", "C", "1")


def make(claim=Claim.CAUSAL, prior=0.0):
    return Hypothesis(
        hid="H", statement_nl="t",
        formal=Triple(R.resolve("GATA6", EntityType.GENE), "regulates",
                      R.resolve("PAAD", EntityType.DISEASE), Direction.UP),
        scope=Scope(disease=R.resolve("PAAD", EntityType.DISEASE)),
        claim=claim, log_prior_odds=prior,
        predictions=[Prediction("p", "d", "m", ">=", 1.0, "assay")],
    )


def ev(kind, lr=50.0, cls="c", **kw):
    return EvidenceItem(eid=f"e{id(kw)}", kind=kind, statement="s", raw_lr=lr,
                        source_id="src", independence_class=cls, **kw)


class TestAttenuation(unittest.TestCase):
    def test_ceiling_is_binding(self):
        for r in (0.5, 0.75, 0.9, 0.95):
            self.assertAlmostEqual(attenuate_lr(1e9, r), information_ceiling(r), places=3)
            self.assertLess(attenuate_lr(1e9, r), information_ceiling(r) + 1e-6)

    def test_identity_and_limits(self):
        self.assertAlmostEqual(attenuate_lr(7.0, 1.0), 7.0)
        self.assertAlmostEqual(attenuate_lr(7.0, 0.0), 1.0)
        self.assertAlmostEqual(attenuate_lr(1.0, 0.6), 1.0)

    def test_symmetric_floor(self):
        r = 0.9
        self.assertAlmostEqual(attenuate_lr(1e-9, r), 1.0 - r, places=6)


class TestDependence(unittest.TestCase):
    def test_redundant_items_do_not_sum(self):
        dep = DependenceModel()
        same = [ev(EvidenceKind.LITERATURE, cls="X") for _ in range(10)]
        diff = [ev(EvidenceKind.LITERATURE, cls=f"X{i}") for i in range(10)]
        s_same, _ = aggregate_log_lr(same, dep)
        s_diff, _ = aggregate_log_lr(diff, dep)
        self.assertLess(s_same, s_diff / 3)


class TestGuards(unittest.TestCase):
    def test_G1_no_prediction_no_registration(self):
        h = make(); h.predictions = []
        with self.assertRaises(GuardError):
            h.register(Ledger(), GEN)

    def test_G2_confirmatory_before_prereg_blocked(self):
        h, L = make(), Ledger()
        with self.assertRaises(GuardError):
            h.attach_evidence(ev(EvidenceKind.OWN_OMICS, confirmatory=True), L, GEN)

    def test_G3_literature_only_never_survives(self):
        h, L = make(Claim.ASSOCIATIVE, prior=2.0), Ledger()
        h.register(L, GEN)
        for i in range(6):
            h.attach_evidence(ev(EvidenceKind.LITERATURE, cls=f"c{i}"), L, GEN)
        for i, a in enumerate((CRIT1, CRIT2, CRIT1)):
            h.record_refutation(RefutationAttempt(f"r{i}", a, "probe", "g", False), L)
        self.assertGreater(h.posterior, 0.9)
        self.assertIs(h.adjudicate(L, CRIT1), State.UNDER_TEST)

    def test_G4_causal_requires_intervention(self):
        h, L = make(Claim.CAUSAL, prior=2.0), Ledger()
        h.register(L, GEN)
        for i in range(4):
            h.attach_evidence(ev(EvidenceKind.OWN_OMICS, cls=f"c{i}"), L, GEN)
        for i, a in enumerate((CRIT1, CRIT2)):
            h.record_refutation(RefutationAttempt(f"r{i}", a, "probe", "g", False), L)
        self.assertGreater(h.posterior, 0.9)
        self.assertIs(h.adjudicate(L, CRIT1), State.UNDER_TEST)
        h.attach_evidence(ev(EvidenceKind.PERTURBATION, cls="pert"), L, GEN)
        self.assertIs(h.adjudicate(L, CRIT1), State.SURVIVING)

    def test_G5_needs_independent_critics(self):
        h, L = make(Claim.ASSOCIATIVE, prior=2.0), Ledger()
        h.register(L, GEN)
        h.attach_evidence(ev(EvidenceKind.PERTURBATION, cls="p"), L, GEN)
        h.record_refutation(RefutationAttempt("r0", CRIT1, "probe", "g", False), L)
        self.assertIs(h.adjudicate(L, CRIT1), State.UNDER_TEST)   # 1 crítico só
        h.record_refutation(RefutationAttempt("r1", CRIT2, "probe", "g", False), L)
        self.assertIs(h.adjudicate(L, CRIT1), State.SURVIVING)

    def test_successful_refutation_kills(self):
        h, L = make(Claim.ASSOCIATIVE, prior=3.0), Ledger()
        h.register(L, GEN)
        h.attach_evidence(ev(EvidenceKind.PERTURBATION, cls="p"), L, GEN)
        h.record_refutation(RefutationAttempt("r", CRIT1, "negative_control", "g", True), L)
        self.assertIs(h.adjudicate(L, CRIT1), State.REFUTED)

    def test_decisive_evidence_kills_regardless_of_posterior(self):
        h, L = make(Claim.ASSOCIATIVE, prior=5.0), Ledger()
        h.register(L, GEN)
        h.attach_evidence(ev(EvidenceKind.WETLAB, lr=1.0, cls="w", decisive=True), L, GEN)
        self.assertIs(h.adjudicate(L, CRIT1), State.REFUTED)


class TestLedger(unittest.TestCase):
    def test_tampering_detected(self):
        h, L = make(), Ledger()
        h.register(L, GEN)
        h.attach_evidence(ev(EvidenceKind.OWN_OMICS, confirmatory=True), L, GEN)
        self.assertTrue(L.verify()[0])
        list(L)[1].payload["raw_lr"] = 9999.0      # adulteração
        self.assertFalse(L.verify()[0])


class TestScheduler(unittest.TestCase):
    def test_eig_peaks_at_uncertainty(self):
        a = Assay("a", "l", EvidenceKind.PERTURBATION, 0.8, 0.85, 1000, 1)
        self.assertGreater(expected_information_gain(0.5, a),
                           expected_information_gain(0.98, a))
        self.assertGreater(expected_information_gain(0.5, a),
                           expected_information_gain(0.02, a))

    def test_uninformative_assay_gives_zero(self):
        useless = Assay("u", "l", EvidenceKind.SIMULATION, 0.5, 0.5, 1, 1)
        self.assertAlmostEqual(expected_information_gain(0.4, useless), 0.0, places=9)


if __name__ == "__main__":
    unittest.main(verbosity=2)
