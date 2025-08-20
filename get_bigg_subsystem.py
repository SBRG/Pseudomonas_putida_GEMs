import requests
from difflib import SequenceMatcher
from collections import Counter
from typing import List, Dict, Tuple, Optional

BIGG_API = "https://bigg.ucsd.edu/api/v2"

# --- HTTP helpers ------------------------------------------------------------

def _get(url, **kwargs):
    """GET with clearer errors and small retry on 5xx."""
    last = None
    for _ in range(2):
        r = requests.get(url, timeout=30, **kwargs)
        if r.status_code >= 500:
            last = r
            continue
        if r.status_code == 429:
            last = r
            continue
        r.raise_for_status()
        return r
    if last is not None:
        raise requests.HTTPError(f"{last.status_code} for {last.url}")
    raise RuntimeError("Unknown HTTP error")

# --- Universal reaction id resolution (no /search) ---------------------------

def _exact_universal_rxn(rxn_id: str) -> Optional[str]:
    """Return rxn_id if it exists as a universal reaction, else None."""
    url = f"{BIGG_API}/universal/reactions/{rxn_id}"
    r = requests.get(url, timeout=20)
    if r.status_code == 200:
        j = r.json()
        return j.get("bigg_id") or rxn_id
    return None

def _all_universal_rxn_ids(max_pages: int = 1000) -> List[str]:
    """
    Page through /universal/reactions to get all reaction IDs.
    Avoids /search entirely.
    """
    ids = []
    page = 1
    while page <= max_pages:
        url = f"{BIGG_API}/universal/reactions"
        r = _get(url, params={"page": page})
        data = r.json()
        items = data.get("results") or data.get("reactions") or []
        # BiGG v2 returns a list of dicts with 'bigg_id'
        page_ids = [it.get("bigg_id") for it in items if it.get("bigg_id")]
        if not page_ids:
            break
        ids.extend(page_ids)
        # stop if we've reached the last page (BiGG doesn't always expose 'next')
        if len(page_ids) < 1000:  # heuristic; BiGG commonly uses page size 1000
            break
        page += 1
    return ids

def _best_universal_reaction_match(query: str, min_sim: float = 0.6) -> Tuple[Optional[str], float]:
    """
    Try exact universal first; else fuzzy-match against the full universal list.
    """
    exact = _exact_universal_rxn(query)
    if exact:
        return exact, 1.0

    # Build a local catalog and fuzzy match
    catalog = _all_universal_rxn_ids()
    if not catalog:
        return None, 0.0

    best = None
    best_score = 0.0
    q = query.lower()
    for bid in catalog:
        s = SequenceMatcher(None, q, bid.lower()).ratio()
        if s > best_score:
            best, best_score = bid, s
    if best_score < min_sim:
        return None, best_score
    return best, best_score

# --- Model enumeration & subsystem lookup (no /search) -----------------------

# A small, biased shortlist to try first (fast) — tweak to your organism(s)
PREFERRED_MODELS = [
    # E. coli (common references)
    "iML1515", "iJO1366", "iAF1260",
    # Pseudomonas putida (frequently used)
    "iJN1462", "iJN1411", "iJP815",
    # Generic/other popular
    "iCHOv1", "Recon3D"
]

def _fetch_model_list(max_models: int = 500) -> List[str]:
    """
    Pull list of model IDs from /models (no paging param documented; do one hit).
    """
    url = f"{BIGG_API}/models"
    r = _get(url)
    data = r.json()
    models = [m.get("bigg_id") for m in data.get("results", []) if m.get("bigg_id")]
    # Keep a reasonable cap
    return models[:max_models]

def _fetch_model_reaction_detail(model_id: str, rxn_id: str) -> Optional[Dict]:
    """
    GET /models/{model_id}/reactions/{rxn_id}
    Returns JSON dict or None on 404/non-200.
    """
    url = f"{BIGG_API}/models/{model_id}/reactions/{rxn_id}"
    r = requests.get(url, timeout=20)
    if r.status_code != 200:
        return None
    return r.json()

def _collect_subsystems_for_rxn(
    rxn_id: str,
    preferred_first: bool = True,
    max_models_probe: int = 60,
) -> Tuple[List[Dict], Counter]:
    """
    Try PREFERRED_MODELS first, then broaden to the global model list until we
    gather enough hits (capped by max_models_probe).
    Returns (per_model, counts) where:
      per_model = [{'model': ..., 'subsystem': ..., 'rxn_name': ...}, ...]
      counts = Counter of subsystem strings (non-empty)
    """
    per_model = []
    counts = Counter()

    tried = set()
    probes = 0

    def try_models(models: List[str]):
        nonlocal probes
        for m in models:
            if probes >= max_models_probe:
                return
            if m in tried:
                continue
            tried.add(m)
            info = _fetch_model_reaction_detail(m, rxn_id)
            probes += 1
            if not info:
                continue
            subsystem = (info.get("subsystem") or "").strip() or None
            name = info.get("name") or None
            per_model.append({"model": m, "subsystem": subsystem, "rxn_name": name})
            if subsystem:
                counts[subsystem] += 1

    # 1) preferred shortlist
    if preferred_first:
        try_models(PREFERRED_MODELS)

    # 2) broaden to all models (in a stable order but skipping already tried)
    if probes < max_models_probe and counts.total() == 0:
        # only fetch model list if we still need more coverage
        all_models = _fetch_model_list(max_models=500)
        # put untried ones first
        rest = [m for m in all_models if m not in tried]
        try_models(rest)

    return per_model, counts

# --- Public function ---------------------------------------------------------

def get_bigg_subsystem(
    reaction_id: str,
    min_similarity: float = 0.6,
    max_models_to_check: int = 60,
    choose_strategy: str = "majority",
) -> Dict:
    """
    Given a reaction identifier (exact or approximate), resolve to a BiGG universal reaction
    WITHOUT using /search, then gather subsystem annotations across models and
    return a consensus.

    Returns a dict:
    {
      'query': str,
      'matched_bigg_reaction': Optional[str],
      'similarity': float,
      'chosen_subsystem': Optional[str],
      'subsystem_counts': Dict[str, int],
      'per_model': List[{'model','subsystem','rxn_name'}]
    }
    """
    rxn_id, sim = _best_universal_reaction_match(reaction_id, min_sim=min_similarity)
    if not rxn_id:
        return {
            "query": reaction_id,
            "matched_bigg_reaction": None,
            "similarity": sim,
            "chosen_subsystem": None,
            "subsystem_counts": {},
            "per_model": [],
        }

    per_model, counts = _collect_subsystems_for_rxn(rxn_id, max_models_probe=max_models_to_check)

    chosen = None
    if choose_strategy == "first_nonempty":
        for row in per_model:
            if row["subsystem"]:
                chosen = row["subsystem"]
                break
    else:
        if counts:
            max_count = max(counts.values())
            chosen = sorted([k for k, v in counts.items() if v == max_count])[0]

    return {
        "query": reaction_id,
        "matched_bigg_reaction": rxn_id,
        "similarity": sim,
        "chosen_subsystem": chosen,
        "subsystem_counts": dict(counts),
        "per_model": per_model,
    }


