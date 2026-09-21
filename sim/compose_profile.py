#!/usr/bin/env python3
"""Composition of named observation channels into a Disease Profile.

Research toy. Seed 20260922. Not a medical device, not a dose, not a cure.
The legal operator is disjoint union. Identification, a product, a mean,
promotion into theta, self-composition, and a cross-class join are refused.
"""

from __future__ import annotations

import copy
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml
from jsonschema import Draft7Validator

ROOT = Path(__file__).resolve().parents[1]
SIM = Path(__file__).resolve().parent
FIG = SIM / "figures"
PROF = ROOT / "profiles"
SCHEMA_PATH = ROOT / "schema" / "composed_profile.schema.json"

SEED = 20260922
SIGMA = 0.06
N_REP = 4
CHI2_95_1DF = 3.841
FLAT_SPREAD_BELOW = 0.5
RANK_TOL = 1e-8
GRID_LO = 0.40
GRID_HI = 2.50
GRID_N = 21

ALPHA = 0.62
BETA = 1.15
GAMMA = 2.40
D_ALPHA = np.array([0.25, 0.50, 0.75, 1.00, 1.25], dtype=float)
D_BETA = np.array([1.80, 0.40, 1.80, 0.40, 1.20], dtype=float)
D_GAMMA = np.array([0.30, 0.90, 0.30, 1.50, 0.60], dtype=float)

DISEASE_ID = "named-channel-class-unspecified"
DISCLAIMER = (
    "Research object. Not a medical device, not clinical decision support, "
    "not a diagnostic, not a dose, and not a cure. No document DOI."
)
SYMBOL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")
SNAPSHOT_DOIS = {
    "10.1093/bioinformatics/btp358",
    "10.1073/pnas.1510507113",
    "10.1038/sdata.2016.18",
    "10.1016/0025-5564(70)90132-x",
    "10.1109/3468.477860",
    "10.1016/j.mbs.2014.08.008",
    "10.1098/rsta.1922.0009",
    "10.1007/978-1-4612-0919-5_16",
}
REQUIRED_NON_PARAMETERS = [
    "merged_therapeutic_theta",
    "composition_quotient",
    "product_of_map_coefficients",
    "arithmetic_mean_of_channels",
    "channel_promoted_into_theta",
    "cross_disease_join",
    "self_composition",
    "assay_scalar_as_treatment",
]
PRIOR_SOURCES = {"T08", "T15", "declared_map"}


def rnd(value):
    if isinstance(value, dict):
        return {k: rnd(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [rnd(v) for v in value]
    if isinstance(value, (np.floating, float)):
        if not np.isfinite(value):
            return None
        return float(f"{float(value):.8g}")
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, np.ndarray):
        return rnd(value.tolist())
    return value


def energy(design: np.ndarray) -> float:
    return float(np.dot(design, design))


def fisher_scale() -> float:
    return N_REP / (SIGMA ** 2)


def legal_information():
    scale = fisher_scale()
    diag = {
        "alpha": scale * energy(D_ALPHA),
        "beta": scale * energy(D_BETA),
        "gamma": scale * energy(D_GAMMA),
    }
    names = ["alpha", "beta", "gamma"]
    matrix = np.diag([diag[name] for name in names])
    return diag, matrix


def rank_of(matrix: np.ndarray) -> int:
    eigenvalues = np.linalg.eigvalsh(np.atleast_2d(matrix))
    eigenvalues = np.real(eigenvalues)
    largest = float(np.max(np.abs(eigenvalues))) if eigenvalues.size else 0.0
    if largest == 0.0:
        return 0
    return int(np.sum(eigenvalues > RANK_TOL * largest))


def spectrum(name: str, matrix: np.ndarray) -> dict:
    eigenvalues = np.sort(np.real(np.linalg.eigvalsh(np.atleast_2d(np.array(matrix, dtype=float)))))[::-1]
    return {
        "name": name,
        "n_param": int(np.atleast_2d(matrix).shape[0]),
        "eigenvalues": eigenvalues,
        "rank": rank_of(matrix),
    }


def identify_estimate(diag: dict) -> float:
    num = ALPHA * diag["alpha"] + BETA * diag["beta"] + GAMMA * diag["gamma"]
    den = diag["alpha"] + diag["beta"] + diag["gamma"]
    return num / den


def product_information() -> np.ndarray:
    coupling = D_ALPHA * D_GAMMA
    sens_alpha = GAMMA * coupling
    sens_gamma = ALPHA * coupling
    scale = fisher_scale()
    matrix = np.zeros((2, 2))
    matrix[0, 0] = scale * float(np.dot(sens_alpha, sens_alpha))
    matrix[1, 1] = scale * float(np.dot(sens_gamma, sens_gamma))
    matrix[0, 1] = matrix[1, 0] = scale * float(np.dot(sens_alpha, sens_gamma))
    return matrix


def relative_grid() -> np.ndarray:
    grid = np.geomspace(GRID_LO, GRID_HI, GRID_N)
    if not np.any(np.isclose(grid, 1.0)):
        grid = np.sort(np.append(grid, 1.0))
    return grid


def call_profile(delta: np.ndarray) -> str:
    spread = float(np.max(delta) - np.min(delta))
    if spread < FLAT_SPREAD_BELOW:
        return "flat"
    if float(delta[0]) > CHI2_95_1DF and float(delta[-1]) > CHI2_95_1DF:
        return "closed"
    return "open"


def profile_record(name: str, grid: np.ndarray, delta: np.ndarray, center: float) -> dict:
    return {
        "name": name,
        "call": call_profile(delta),
        "spread": float(np.max(delta) - np.min(delta)),
        "center": center,
        "endpoint_grid": [float(grid[0]), float(grid[-1])],
        "endpoint_delta_chi2": [float(delta[0]), float(delta[-1])],
        "argmin": float(grid[int(np.argmin(delta))]),
        "min_delta_chi2": float(np.min(delta)),
    }


def non_parameter_catalog() -> list:
    text = {
        "merged_therapeutic_theta": "One symbol used as every channel's map coefficient and entered as a therapeutic parameter.",
        "composition_quotient": "A join that identifies distinct map coefficients.",
        "product_of_map_coefficients": "A product of two map coefficients written as one parameter.",
        "arithmetic_mean_of_channels": "An arithmetic mean of channel readouts written as one parameter.",
        "channel_promoted_into_theta": "A single named channel whose status is changed to in_theta by the act of joining.",
        "cross_disease_join": "A disjoint union across two disease-class identifiers.",
        "self_composition": "A channel joined to a second copy of the same channel identifier.",
        "assay_scalar_as_treatment": "An assay-shaped map coefficient copied in as a treatment parameter. No assay table is stored here.",
    }
    return [{"id": key, "statement": text[key]} for key in REQUIRED_NON_PARAMETERS]


def citations() -> list:
    return [
        {
            "id": "cit-raue-2009",
            "vancouver": "Raue A, Kreutz C, Maiwald T, Bachmann J, Schilling M, Klingmüller U, et al. Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. Bioinformatics. 2009;25(15):1923-1929.",
            "doi": "10.1093/bioinformatics/btp358",
        },
        {
            "id": "cit-bareinboim-2016",
            "vancouver": "Bareinboim E, Pearl J. Causal inference and the data-fusion problem. Proc Natl Acad Sci U S A. 2016;113(27):7345-7352.",
            "doi": "10.1073/pnas.1510507113",
        },
        {
            "id": "cit-wilkinson-2016",
            "vancouver": "Wilkinson MD, Dumontier M, Aalbersberg IJ, Appleton G, Axton M, Baak A, et al. The FAIR Guiding Principles for scientific data management and stewardship. Sci Data. 2016;3:160018.",
            "doi": "10.1038/sdata.2016.18",
        },
        {
            "id": "cit-bellman-1970",
            "vancouver": "Bellman R, Åström KJ. On structural identifiability. Math Biosci. 1970;7(3-4):329-339.",
            "doi": "10.1016/0025-5564(70)90132-x",
        },
        {
            "id": "cit-bloch-1996",
            "vancouver": "Bloch I. Information combination operators for data fusion: a comparative review with classification. IEEE Trans Syst Man Cybern A. 1996;26(1):52-67.",
            "doi": "10.1109/3468.477860",
        },
        {
            "id": "cit-eisenberg-2014",
            "vancouver": "Eisenberg MC, Hayashi MAL. Determining identifiable parameter combinations using subset profiling. Math Biosci. 2014;256:116-126.",
            "doi": "10.1016/j.mbs.2014.08.008",
        },
        {
            "id": "cit-fisher-1922",
            "vancouver": "Fisher RA. On the mathematical foundations of theoretical statistics. Philos Trans R Soc Lond A. 1922;222:309-368.",
            "doi": "10.1098/rsta.1922.0009",
        },
        {
            "id": "cit-rao-1992",
            "vancouver": "Rao CR. Information and the accuracy attainable in the estimation of statistical parameters. In: Breakthroughs in statistics. Springer Series in Statistics. 1992.",
            "doi": "10.1007/978-1-4612-0919-5_16",
        },
        {
            "id": "cit-t03",
            "vancouver": "Ogbonna KE. Disease profiles for complex pathologies: a gated method for systemic personalized-medicine research objects [Internet]. Thesis #3. 2026 [cited 2026 Sep 21]. Available from: https://github.com/cloudynirvana/thesis-03-disease-profile",
        },
        {
            "id": "cit-t08",
            "vancouver": "Ogbonna KE. Green-synthesized silver nanoparticles from Carica papaya as an in-vitro metabolic observation channel [Internet]. Thesis #8. 2026 [cited 2026 Sep 21]. Available from: https://github.com/cloudynirvana/thesis-08-papaya-agnp-observation-channel",
        },
        {
            "id": "cit-t15",
            "vancouver": "Ogbonna KE. Encoding a compartmental AgNP-exosome-Raman theranostic concept as a multi-observation Disease Profile research object [Internet]. Thesis #15. 2026 [cited 2026 Sep 21]. Available from: https://github.com/cloudynirvana/thesis-15-nanobiocomposite-multiobservation-profile",
        },
    ]


def channel(channel_id, name, symbol, coefficient, observes, factor_id, source, refusal) -> dict:
    return {
        "channel_id": channel_id,
        "name": name,
        "symbol": symbol,
        "map_coefficient": coefficient,
        "observes": observes,
        "parameter_status": "not_in_theta",
        "factor_id": factor_id,
        "source_prior": source,
        "refusal": refusal,
    }


def amyl_channel(factor_id="f_amyl") -> dict:
    return channel(
        "ch_amyl_map",
        "Amylase-shaped observation map",
        "y_A",
        "alpha",
        "A declared saturating function of a substrate coordinate. The coefficient is not an inhibition percentage and not an IC50.",
        factor_id,
        "T08",
        "The map coefficient stays out of theta. Joining this channel to other maps does not convert it into a treatment parameter.",
    )


def spectral_channel(factor_id="f_spectral") -> dict:
    return channel(
        "ch_spectral",
        "Spectral-contrast map",
        "y_S",
        "beta",
        "A declared spectral contrast. Not a copy of a published reporter loading.",
        factor_id,
        "T15",
        "The spectral coefficient is not a killing coefficient and is not renamed when a second channel is attached.",
    )


def compartment_channel(factor_id="f_compartment") -> dict:
    return channel(
        "ch_compartment",
        "Compartment-indicator map",
        "y_C",
        "gamma",
        "A declared compartment indicator. Not a vesicle count and not a delivery rate.",
        factor_id,
        "T15",
        "The compartment coefficient stays a map coefficient. Tropism language does not become a rate inside theta.",
    )


def answers() -> list:
    return [
        {
            "id": "Q1",
            "question": "Who is asking?",
            "label": "A methods researcher asking when named observation channels may be joined.",
        },
        {
            "id": "Q2",
            "question": "Which disease, not a generic slogan?",
            "label": "An unspecified class of named observation channels. No person and no stage.",
        },
        {
            "id": "Q3",
            "question": "What is the current regime?",
            "label": "An in-silico join of declared maps. No clinical regime is asserted.",
        },
        {
            "id": "Q4",
            "question": "Where is the system stuck?",
            "label": "The portfolio has a profile object, one named channel, and one multi-channel file. It does not yet say which joins preserve that file as a profile.",
        },
    ]


def observables() -> list:
    return [
        {
            "id": "obs-map-not-parameter",
            "statement": "A named observation channel is a map. Its coefficient is not a member of a therapeutic parameter list unless a separate identification argument says so.",
            "citation_ids": ["cit-bellman-1970", "cit-t08"],
        },
        {
            "id": "obs-profile-object",
            "statement": "A Disease Profile is a versioned research export. It is not a chart and not a parameter table.",
            "citation_ids": ["cit-t03", "cit-wilkinson-2016"],
        },
        {
            "id": "obs-fusion-is-not-join",
            "statement": "Pooling readings under one symbol answers a different question from keeping the readings as separate maps.",
            "citation_ids": ["cit-bareinboim-2016", "cit-bloch-1996", "cit-t15"],
        },
    ]


def mechanisms(merged_status="refused", falsifier_merged="") -> list:
    return [
        {
            "id": "M-DISJOINT",
            "statement": "The composed object keeps one symbol and one map coefficient for each named channel.",
            "evidence_class": "definition",
            "falsifier": "A preparation series whose loadings differ, fitted anyway by one shared coefficient.",
            "status": "admitted",
        },
        {
            "id": "M-ONE-THETA",
            "statement": "One therapeutic parameter accounts for every channel in the join.",
            "evidence_class": "knowledge",
            "falsifier": falsifier_merged,
            "status": merged_status,
            "refusal_rule": "composition_quotient" if merged_status == "refused" else "",
        },
    ]


def hypotheses(status="forbidden_to_enter_theta", statement=None) -> list:
    if statement is None:
        statement = (
            "A disjoint union of named maps is a Disease Profile only while each map "
            "coefficient stays out of theta and the composition record names no quotient."
        )
    return [
        {
            "id": "H-JOIN-001",
            "statement": statement,
            "falsifier": "Any composed file with a non-empty theta list, a quotient, a combiner, a repeated channel identifier, or two disease-class identifiers.",
            "parameter_status": status,
        }
    ]


def composition(operator, factors, quotient=None, combiner=None) -> dict:
    return {
        "operator": operator,
        "quotient": quotient,
        "combiner": combiner,
        "factors": factors,
    }


def factor(factor_id, source, disease_id, channel_ids) -> dict:
    return {
        "factor_id": factor_id,
        "source_prior": source,
        "disease_id": disease_id,
        "channel_ids": list(channel_ids),
    }


def profile_shell(profile_id, channels, comp, disease_id=DISEASE_ID, disease_label=None) -> dict:
    if disease_label is None:
        disease_label = (
            "Unspecified class used to test composition of named maps. Not a patient and not a stage."
        )
    return {
        "schema_version": "1.2.0-compose",
        "profile_id": profile_id,
        "extends": {
            "object": "DiseaseProfile",
            "schema_version": "1.0.0",
            "source": "Thesis #3 Disease Profile contract, local composition extension 1.2.0-compose. Does not edit schema 1.1.0-multiobs.",
        },
        "disease_id": disease_id,
        "disease_label": disease_label,
        "created_at": "2026-09-21T00:00:00Z",
        "asker_role": "methods_researcher",
        "answers": answers(),
        "observables": observables(),
        "observation_channels": channels,
        "composition": comp,
        "candidate_mechanisms": mechanisms(),
        "non_parameters": non_parameter_catalog(),
        "theta": [],
        "admitted_hypotheses": hypotheses(),
        "citations": citations(),
        "disclaimer": DISCLAIMER,
    }


def legal_bundle() -> dict:
    channels = [
        amyl_channel("f_amyl"),
        spectral_channel("f_t15_pair"),
        compartment_channel("f_t15_pair"),
    ]
    factors = [
        factor("f_amyl", "T08", DISEASE_ID, ["ch_amyl_map"]),
        factor("f_t15_pair", "T15", DISEASE_ID, ["ch_spectral", "ch_compartment"]),
    ]
    return profile_shell(
        "dp-compose-bundle-001",
        channels,
        composition("disjoint_union", factors),
    )


def legal_singletons() -> dict:
    channels = [
        amyl_channel("f_amyl"),
        spectral_channel("f_spectral"),
        compartment_channel("f_compartment"),
    ]
    factors = [
        factor("f_amyl", "T08", DISEASE_ID, ["ch_amyl_map"]),
        factor("f_spectral", "T15", DISEASE_ID, ["ch_spectral"]),
        factor("f_compartment", "T15", DISEASE_ID, ["ch_compartment"]),
    ]
    doc = profile_shell(
        "dp-compose-singletons-001",
        channels,
        composition("disjoint_union", factors),
    )
    return doc


def legal_binary() -> dict:
    channels = [
        spectral_channel("f_spectral"),
        compartment_channel("f_compartment"),
    ]
    factors = [
        factor("f_spectral", "T15", DISEASE_ID, ["ch_spectral"]),
        factor("f_compartment", "T15", DISEASE_ID, ["ch_compartment"]),
    ]
    return profile_shell(
        "dp-compose-binary-001",
        channels,
        composition("disjoint_union", factors),
    )


def illegal_identify() -> dict:
    doc = legal_singletons()
    doc["profile_id"] = "dp-compose-illegal-identify-001"
    for item in doc["observation_channels"]:
        item["symbol"] = "theta_theranostic"
        item["map_coefficient"] = "theta_theranostic"
        item["parameter_status"] = "in_theta"
        item["refusal"] = "Removed. The channel is the parameter."
    doc["composition"]["operator"] = "identify"
    doc["composition"]["quotient"] = {
        "identifies": ["alpha", "beta", "gamma"],
        "as": "theta_theranostic",
    }
    doc["composition"]["combiner"] = "identify"
    doc["theta"] = ["theta_theranostic"]
    doc["candidate_mechanisms"] = mechanisms("admitted", "")
    doc["admitted_hypotheses"] = hypotheses(
        "entered_theta",
        "A single theranostic parameter accounts for every named map in the join.",
    )
    doc["non_parameters"] = [
        item for item in doc["non_parameters"] if item["id"] != "merged_therapeutic_theta"
    ]
    return doc


def illegal_product() -> dict:
    doc = legal_singletons()
    doc["profile_id"] = "dp-compose-illegal-product-001"
    doc["composition"]["operator"] = "quotient"
    doc["composition"]["quotient"] = {"of": ["alpha", "gamma"], "as": "pi_alpha_gamma"}
    doc["composition"]["combiner"] = "product"
    doc["theta"] = ["pi_alpha_gamma"]
    doc["candidate_mechanisms"] = mechanisms("admitted", "")
    doc["admitted_hypotheses"] = hypotheses(
        "entered_theta",
        "The product of the amylase-shaped coefficient and the compartment coefficient is the therapeutic parameter.",
    )
    doc["non_parameters"] = [
        item for item in doc["non_parameters"] if item["id"] != "product_of_map_coefficients"
    ]
    return doc


def illegal_mean() -> dict:
    doc = legal_singletons()
    doc["profile_id"] = "dp-compose-illegal-mean-001"
    doc["composition"]["operator"] = "arithmetic_mean"
    doc["composition"]["combiner"] = "arithmetic_mean"
    doc["composition"]["quotient"] = {"of": ["alpha", "beta", "gamma"], "as": "theta_mean"}
    doc["theta"] = ["theta_mean"]
    doc["candidate_mechanisms"] = mechanisms("admitted", "")
    doc["admitted_hypotheses"] = hypotheses(
        "entered_theta",
        "The arithmetic mean of the three map coefficients is the therapeutic parameter.",
    )
    doc["non_parameters"] = [
        item for item in doc["non_parameters"] if item["id"] != "arithmetic_mean_of_channels"
    ]
    return doc


def illegal_promote() -> dict:
    doc = legal_bundle()
    doc["profile_id"] = "dp-compose-illegal-promote-001"
    for item in doc["observation_channels"]:
        if item["channel_id"] == "ch_amyl_map":
            item["parameter_status"] = "in_theta"
            item["refusal"] = "Removed. Joining promoted this map into theta."
    doc["theta"] = ["alpha"]
    doc["candidate_mechanisms"] = mechanisms("admitted", "")
    doc["admitted_hypotheses"] = hypotheses(
        "entered_theta",
        "Attaching the amylase-shaped map to the other channels writes alpha into theta.",
    )
    doc["non_parameters"] = [
        item for item in doc["non_parameters"] if item["id"] != "channel_promoted_into_theta"
    ]
    return doc


def illegal_self() -> dict:
    doc = legal_binary()
    doc["profile_id"] = "dp-compose-illegal-self-001"
    duplicate = copy.deepcopy(doc["observation_channels"][0])
    duplicate["name"] = "Second copy of the spectral map"
    doc["observation_channels"].append(duplicate)
    doc["composition"]["factors"].append(
        factor("f_spectral_again", "T15", DISEASE_ID, ["ch_spectral"])
    )
    doc["non_parameters"] = [
        item for item in doc["non_parameters"] if item["id"] != "self_composition"
    ]
    return doc


def illegal_cross() -> dict:
    doc = legal_binary()
    doc["profile_id"] = "dp-compose-illegal-cross-001"
    doc["composition"]["factors"][1]["disease_id"] = "other-channel-class"
    doc["non_parameters"] = [
        item for item in doc["non_parameters"] if item["id"] != "cross_disease_join"
    ]
    return doc


def failed_rules(doc: dict) -> list:
    failed = []
    if doc.get("schema_version") != "1.2.0-compose":
        failed.append("schema_version")
    extends = doc.get("extends") or {}
    if extends.get("object") != "DiseaseProfile" or extends.get("schema_version") != "1.0.0":
        failed.append("extends_parent")
    answer_ids = [item.get("id") for item in doc.get("answers") or []]
    if answer_ids[:4] != ["Q1", "Q2", "Q3", "Q4"]:
        failed.append("four_answers")
    channels = doc.get("observation_channels") or []
    if len(channels) < 2:
        failed.append("channel_count")
    ids = [item.get("channel_id") for item in channels]
    symbols = [item.get("symbol") for item in channels]
    coeffs = [item.get("map_coefficient") for item in channels]
    if len(ids) != len(set(ids)):
        failed.append("channel_id_unique")
    if len(symbols) != len(set(symbols)):
        failed.append("symbol_unique")
    if len(coeffs) != len(set(coeffs)):
        failed.append("map_coefficient_unique")
    if any(not SYMBOL_RE.match(str(item.get("symbol", ""))) for item in channels):
        failed.append("symbol_not_name")
    if any(not SYMBOL_RE.match(str(item.get("map_coefficient", ""))) for item in channels):
        failed.append("map_coefficient_not_name")
    if any(item.get("parameter_status") != "not_in_theta" for item in channels):
        failed.append("channel_parameter_status")
    theta = doc.get("theta")
    if theta != []:
        failed.append("theta_not_empty")
    comp = doc.get("composition") or {}
    if comp.get("operator") != "disjoint_union":
        failed.append("operator_not_disjoint_union")
    if comp.get("quotient") is not None:
        failed.append("quotient_not_null")
    if comp.get("combiner") is not None:
        failed.append("combiner_not_null")
    factors = comp.get("factors") or []
    if len(factors) < 2:
        failed.append("factor_count")
    covered = []
    factor_ids = []
    disease_ids = []
    for item in factors:
        factor_ids.append(item.get("factor_id"))
        disease_ids.append(item.get("disease_id"))
        if item.get("source_prior") not in PRIOR_SOURCES:
            failed.append("factor_source")
        covered.extend(item.get("channel_ids") or [])
    if len(factor_ids) != len(set(factor_ids)):
        failed.append("factor_id_unique")
    if covered != ids or len(covered) != len(set(covered)):
        failed.append("factor_partition")
    profile_disease = doc.get("disease_id")
    if any(item != profile_disease for item in disease_ids) or len(set(disease_ids)) != 1:
        failed.append("disease_mismatch")
    label = str(profile_disease or "").lower()
    if any(token in label for token in ("patient", "mrn", "subject")):
        failed.append("disease_id_names_person")
    present = {item.get("id") for item in doc.get("non_parameters") or []}
    for key in REQUIRED_NON_PARAMETERS:
        if key not in present:
            failed.append(f"missing_non_parameter:{key}")
    merged = next((item for item in doc.get("candidate_mechanisms") or [] if item.get("id") == "M-ONE-THETA"), None)
    if merged is None or merged.get("status") != "refused":
        failed.append("merged_claim_not_refused")
    for item in doc.get("candidate_mechanisms") or []:
        if item.get("status") == "admitted" and not str(item.get("falsifier") or "").strip():
            failed.append("candidate_admitted_without_falsifier")
            break
    for item in doc.get("admitted_hypotheses") or []:
        if item.get("parameter_status") != "forbidden_to_enter_theta":
            failed.append("hypothesis_parameter_status")
            break
        if not str(item.get("falsifier") or "").strip():
            failed.append("hypothesis_missing_falsifier")
            break
    for item in doc.get("citations") or []:
        if "doi" in item:
            doi = str(item.get("doi") or "")
            if not DOI_RE.match(doi):
                failed.append("doi_pattern")
                break
            if doi.lower() not in SNAPSHOT_DOIS:
                failed.append("doi_not_in_snapshot")
                break
    if doc.get("disclaimer") != DISCLAIMER:
        failed.append("disclaimer")
    # Stable unique order.
    seen = set()
    ordered = []
    for item in failed:
        if item not in seen:
            seen.add(item)
            ordered.append(item)
    return ordered


def schema_ok(doc: dict, validator: Draft7Validator) -> bool:
    return not list(validator.iter_errors(doc))


def dump_yaml(path: Path, doc: dict) -> None:
    path.write_text(
        yaml.safe_dump(doc, sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )


def cramer_rao(diag: dict) -> dict:
    out = {}
    truth = {"alpha": ALPHA, "beta": BETA, "gamma": GAMMA}
    for name, info in diag.items():
        se = info ** -0.5
        out[name] = {"se": se, "relative_se": se / truth[name]}
    return out


def delta_at(info: float, grid: np.ndarray, center: float) -> np.ndarray:
    return info * (grid - center) ** 2


def noisy_identify(diag: dict, rng: np.random.Generator) -> dict:
    designs = {"alpha": D_ALPHA, "beta": D_BETA, "gamma": D_GAMMA}
    truth = {"alpha": ALPHA, "beta": BETA, "gamma": GAMMA}
    numer = 0.0
    denom = 0.0
    for name, design in designs.items():
        draws = truth[name] * design[:, None] + rng.normal(0.0, SIGMA, size=(design.size, N_REP))
        numer += float(np.sum(draws * design[:, None]))
        denom += N_REP * energy(design)
    theta_hat = numer / denom
    info = sum(diag.values())
    return {"theta_hat": theta_hat, "information": info}


def plot_spectra(path: Path, spectra: list) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    width = 0.18
    anchors = np.arange(3)
    series = spectra[:7]
    offsets = np.linspace(-0.36, 0.36, len(series))
    for offset, item in zip(offsets, series):
        values = list(item["eigenvalues"])
        while len(values) < 3:
            values.append(0.0)
        ax.bar(anchors + offset, values[:3], width=width, label=item["name"])
    ax.set_xticks(anchors)
    ax.set_xticklabels(["λ1", "λ2", "λ3"])
    ax.set_ylabel("Fisher eigenvalue")
    ax.set_title("Information under join and under quotient")
    ax.legend(fontsize=7, ncol=2, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_profiles(path: Path, curves: list) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    for item in curves:
        ax.plot(item["grid"], item["delta"], label=item["name"])
    ax.axhline(CHI2_95_1DF, color="black", lw=0.8, ls="--", label="χ² 3.841")
    ax.set_xlabel("Coefficient on its own grid")
    ax.set_ylabel("Δχ²")
    ax.set_title("Profiles of a legal coefficient, a product, and a merged symbol")
    ax.legend(fontsize=7, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_ranks(path: Path, rows: list) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    names = [row["name"] for row in rows]
    ranks = [row["rank"] for row in rows]
    colors = ["#1f4e79" if row["accepted"] else "#8c2f39" for row in rows]
    ax.barh(names[::-1], ranks[::-1], color=colors[::-1])
    ax.set_xlabel("Numerical rank of the declared parameter")
    ax.set_title("Accepted joins add rank. Refused quotients do not.")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    PROF.mkdir(parents=True, exist_ok=True)
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    validator = Draft7Validator(schema)

    documents = {
        "legal_bundle": legal_bundle(),
        "legal_singletons": legal_singletons(),
        "legal_binary": legal_binary(),
        "illegal_identify": illegal_identify(),
        "illegal_product": illegal_product(),
        "illegal_mean": illegal_mean(),
        "illegal_promote": illegal_promote(),
        "illegal_self": illegal_self(),
        "illegal_cross": illegal_cross(),
    }
    filenames = {
        "legal_bundle": "legal_bundle.profile.yaml",
        "legal_singletons": "legal_singletons.profile.yaml",
        "legal_binary": "legal_binary.profile.yaml",
        "illegal_identify": "illegal_identify.profile.yaml",
        "illegal_product": "illegal_product.profile.yaml",
        "illegal_mean": "illegal_mean.profile.yaml",
        "illegal_promote": "illegal_promote.profile.yaml",
        "illegal_self": "illegal_self.profile.yaml",
        "illegal_cross": "illegal_cross.profile.yaml",
    }
    reports = {}
    for name, doc in documents.items():
        rules = failed_rules(doc)
        reports[name] = {
            "rules_failed": rules,
            "accepted": rules == [],
            "schema_valid": schema_ok(doc, validator),
        }
        dump_yaml(PROF / filenames[name], doc)
    (PROF / "legal_bundle.profile.json").write_text(
        json.dumps(documents["legal_bundle"], indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    fake = copy.deepcopy(documents["legal_bundle"])
    fake["citations"][0]["doi"] = "10.1234/not-a-real-record"
    fake_failed = failed_rules(fake)

    diag, joint = legal_information()
    amyl_only = np.diag([diag["alpha"], 0.0, 0.0])
    spectral_only = np.diag([0.0, diag["beta"], 0.0])
    compartment_only = np.diag([0.0, 0.0, diag["gamma"]])
    binary = np.diag([0.0, diag["beta"], diag["gamma"]])
    identified = np.array([[sum(diag.values())]])
    product = product_information()
    spectra = [
        spectrum("legal_joint", joint),
        spectrum("amyl_only", amyl_only),
        spectrum("spectral_only", spectral_only),
        spectrum("compartment_only", compartment_only),
        spectrum("binary_spectral_compartment", binary),
        spectrum("illegal_identify", identified),
        spectrum("product_alpha_gamma", product),
    ]

    theta_hat = identify_estimate(diag)
    weights = {name: diag[name] / sum(diag.values()) for name in diag}
    arithmetic_mean = (ALPHA + BETA + GAMMA) / 3.0
    product_value = ALPHA * GAMMA
    product_eigs = np.linalg.eigvalsh(product)
    kernel = np.array([ALPHA, -GAMMA], dtype=float)
    kernel = kernel / np.linalg.norm(kernel)
    # The smallest eigenvector of a rank-1 Gram matrix of parallel columns.
    _, evecs = np.linalg.eigh(product)
    null_vec = evecs[:, 0]
    # Align sign with (alpha, -gamma).
    if np.dot(null_vec, kernel) < 0:
        null_vec = -null_vec

    grid_rel = relative_grid()
    grid_alpha = grid_rel * ALPHA
    grid_beta = grid_rel * BETA
    grid_gamma = grid_rel * GAMMA
    grid_theta = grid_rel * theta_hat
    delta_alpha = delta_at(diag["alpha"], grid_alpha, ALPHA)
    delta_beta = delta_at(diag["beta"], grid_beta, BETA)
    delta_gamma = delta_at(diag["gamma"], grid_gamma, GAMMA)
    delta_theta = delta_at(sum(diag.values()), grid_theta, theta_hat)
    delta_foreign = np.zeros_like(grid_gamma)
    delta_product = np.zeros_like(grid_alpha)
    slice_info = product[0, 0] / (GAMMA ** 2) * (GAMMA ** 2)
    # slice_info computed directly:
    coupling = D_ALPHA * D_GAMMA
    slice_info = fisher_scale() * (GAMMA ** 2) * float(np.dot(coupling, coupling))
    delta_slice = delta_at(slice_info, grid_alpha, ALPHA)

    profiles = {
        "alpha": profile_record("alpha", grid_alpha, delta_alpha, ALPHA),
        "beta": profile_record("beta", grid_beta, delta_beta, BETA),
        "gamma": profile_record("gamma", grid_gamma, delta_gamma, GAMMA),
        "spectral_only_alpha": profile_record("spectral_only_alpha", grid_alpha, np.zeros_like(grid_alpha), ALPHA),
        "amyl_only_gamma": profile_record("amyl_only_gamma", grid_gamma, delta_foreign, GAMMA),
        "product_alpha": profile_record("product_alpha", grid_alpha, delta_product, ALPHA),
        "product_slice_alpha": profile_record("product_slice_alpha", grid_alpha, delta_slice, ALPHA),
        "illegal_theta": profile_record("illegal_theta", grid_theta, delta_theta, theta_hat),
    }
    truth_delta = {
        name: float(sum(diag.values()) * (value - theta_hat) ** 2)
        for name, value in {"alpha": ALPHA, "beta": BETA, "gamma": GAMMA}.items()
    }

    rng = np.random.default_rng(SEED)
    noisy = noisy_identify(diag, rng)
    # One noisy legal coefficient, alpha, same seed stream continued.
    draws = ALPHA * D_ALPHA[:, None] + rng.normal(0.0, SIGMA, size=(D_ALPHA.size, N_REP))
    alpha_hat = float(np.sum(draws * D_ALPHA[:, None]) / (N_REP * energy(D_ALPHA)))
    delta_alpha_noisy = delta_at(diag["alpha"], grid_alpha, alpha_hat)
    profiles["alpha_noisy"] = profile_record("alpha_noisy", grid_alpha, delta_alpha_noisy, alpha_hat)

    bundle_ids = [item["channel_id"] for item in documents["legal_bundle"]["observation_channels"]]
    single_ids = [item["channel_id"] for item in documents["legal_singletons"]["observation_channels"]]
    binary_ids = [item["channel_id"] for item in documents["legal_binary"]["observation_channels"]]
    left = set(["ch_amyl_map"]) | set(["ch_spectral", "ch_compartment"])
    right = set(["ch_amyl_map", "ch_spectral"]) | set(["ch_compartment"])

    algebra = {
        "bundle_equals_singletons": bundle_ids == single_ids,
        "associative_channel_sets": left == right == set(single_ids),
        "commutative_binary": set(binary_ids) == set(reversed(binary_ids)),
        "binary_is_subset": set(binary_ids).issubset(set(single_ids)),
        "rank_adds": spectra[0]["rank"] == spectra[1]["rank"] + spectra[2]["rank"] + spectra[3]["rank"],
        "binary_rank_adds": spectra[4]["rank"] == spectra[2]["rank"] + spectra[3]["rank"],
        "identify_rank": spectra[5]["rank"],
        "product_rank": spectra[6]["rank"],
        "product_null_cosine": float(abs(np.dot(null_vec, kernel))),
        "off_diagonal_max": float(np.max(np.abs(joint - np.diag(np.diag(joint))))),
    }

    rank_rows = [
        {"name": "amyl only", "rank": spectra[1]["rank"], "accepted": True},
        {"name": "spectral only", "rank": spectra[2]["rank"], "accepted": True},
        {"name": "compartment only", "rank": spectra[3]["rank"], "accepted": True},
        {"name": "binary join", "rank": spectra[4]["rank"], "accepted": reports["legal_binary"]["accepted"]},
        {"name": "triple join", "rank": spectra[0]["rank"], "accepted": reports["legal_bundle"]["accepted"]},
        {"name": "identify", "rank": spectra[5]["rank"], "accepted": reports["illegal_identify"]["accepted"]},
        {"name": "product", "rank": spectra[6]["rank"], "accepted": reports["illegal_product"]["accepted"]},
    ]

    plot_spectra(FIG / "fisher_spectra.png", spectra)
    plot_profiles(
        FIG / "profiles_join_quotient.png",
        [
            {"name": "legal α", "grid": grid_rel, "delta": delta_alpha},
            {"name": "product α (hyperbola)", "grid": grid_rel, "delta": delta_product},
            {"name": "product slice, γ fixed", "grid": grid_rel, "delta": delta_slice},
            {"name": "merged θ", "grid": grid_rel, "delta": delta_theta},
        ],
    )
    plot_ranks(FIG / "composition_ranks.png", rank_rows)

    legal_names = ["legal_bundle", "legal_singletons", "legal_binary"]
    illegal_names = [
        "illegal_identify",
        "illegal_product",
        "illegal_mean",
        "illegal_promote",
        "illegal_self",
        "illegal_cross",
    ]
    schema_caught = ["illegal_identify", "illegal_product", "illegal_mean", "illegal_promote"]
    checks = {
        "legal_accepted": all(reports[name]["accepted"] and reports[name]["schema_valid"] for name in legal_names),
        "illegal_refused": all(not reports[name]["accepted"] for name in illegal_names),
        "schema_refuses_quotient_files": all(not reports[name]["schema_valid"] for name in schema_caught),
        "schema_misses_self_and_cross": (
            reports["illegal_self"]["schema_valid"]
            and reports["illegal_cross"]["schema_valid"]
            and not reports["illegal_self"]["accepted"]
            and not reports["illegal_cross"]["accepted"]
        ),
        "fake_doi_refused": "doi_not_in_snapshot" in fake_failed,
        "algebra": algebra["bundle_equals_singletons"] and algebra["associative_channel_sets"] and algebra["rank_adds"],
        "product_kernel": algebra["product_null_cosine"] > 1.0 - 1e-8,
        "profiles": profiles["alpha"]["call"] == "closed"
        and profiles["product_alpha"]["call"] == "flat"
        and profiles["product_slice_alpha"]["call"] == "closed"
        and profiles["illegal_theta"]["call"] == "closed"
        and profiles["amyl_only_gamma"]["call"] == "flat",
    }
    if not all(checks.values()):
        raise SystemExit(f"composition checks failed: {checks}")

    doi_count = sum(1 for item in documents["legal_bundle"]["citations"] if "doi" in item)
    results = {
        "seed": SEED,
        "sigma": SIGMA,
        "n_rep": N_REP,
        "chi2_95_1df": CHI2_95_1DF,
        "flat_spread_below": FLAT_SPREAD_BELOW,
        "rank_tol_relative": RANK_TOL,
        "truth": {"alpha": ALPHA, "beta": BETA, "gamma": GAMMA},
        "design": {
            "alpha": D_ALPHA.tolist(),
            "beta": D_BETA.tolist(),
            "gamma": D_GAMMA.tolist(),
        },
        "design_energy": {name: energy(design) for name, design in {
            "alpha": D_ALPHA, "beta": D_BETA, "gamma": D_GAMMA
        }.items()},
        "information_diagonal": diag,
        "cramer_rao_legal": cramer_rao(diag),
        "identify": {
            "theta_hat": theta_hat,
            "weights": weights,
            "arithmetic_mean_of_coefficients": arithmetic_mean,
            "abs_error": {
                "alpha": abs(theta_hat - ALPHA),
                "beta": abs(theta_hat - BETA),
                "gamma": abs(theta_hat - GAMMA),
            },
            "delta_chi2_at_truth": truth_delta,
            "information": sum(diag.values()),
        },
        "product": {
            "pi": product_value,
            "eigenvalues": np.sort(np.real(product_eigs))[::-1],
            "rank": rank_of(product),
            "null_cosine_with_alpha_neg_gamma": float(abs(np.dot(null_vec, kernel))),
        },
        "noisy": {
            "identify_theta_hat": noisy["theta_hat"],
            "alpha_hat": alpha_hat,
        },
        "validator": reports,
        "fake_doi_rules_failed": fake_failed,
        "algebra": algebra,
        "spectra": spectra,
        "profiles": profiles,
        "profile_doi_count": doi_count,
        "checks": checks,
    }
    (SIM / "results.json").write_text(
        json.dumps(rnd(results), indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(rnd({
        "checks": checks,
        "identify_theta_hat": theta_hat,
        "weights": weights,
        "information_diagonal": diag,
        "product_pi": product_value,
        "profiles": {key: {"call": value["call"], "endpoint_delta_chi2": value["endpoint_delta_chi2"]} for key, value in profiles.items()},
        "truth_delta": truth_delta,
        "cramer_rao": cramer_rao(diag),
        "validator_failed": {name: reports[name]["rules_failed"] for name in illegal_names},
        "spectra": spectra,
        "noisy": results["noisy"],
        "abs_error": results["identify"]["abs_error"],
        "design_energy": results["design_energy"],
    }), indent=2))


if __name__ == "__main__":
    main()
