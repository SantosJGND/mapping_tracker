import json
import os


def _env_int(name, default):
    value = os.environ.get(name)
    if value is None:
        return default
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if parsed < 1:
        raise ValueError(f"{name} must be >= 1")
    return parsed


def _defaults():
    return {
        "workers": _env_int("MAPPING_WORKERS", 4),
        "threads": _env_int("MAPPING_THREADS", max(1, os.cpu_count() or 1)),
        "dedup_heap": os.environ.get("MAPPING_DEDUP_HEAP", "8g"),
        "sort_memory": os.environ.get("MAPPING_SORT_MEMORY"),
        "no_filter": False,
        "no_prefilter": False,
        "no_dedup": False,
        "no_trim": False,
    }


NEGATIVE_KEYS = ("no_filter", "no_prefilter", "no_dedup", "no_trim")
POSITIVE_TO_NEGATIVE = {
    "filter": "no_filter",
    "prefilter": "no_prefilter",
    "dedup": "no_dedup",
    "trim": "no_trim",
}


def resolved_config_path(explicit=None):
    return explicit or os.environ.get("MAPPING_CONFIG") or "config/mapping.json"


def _validate(config, path):
    for key in ("workers", "threads"):
        value = config[key]
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{path}: {key} must be an integer >= 1")
    for key in ("dedup_heap", "sort_memory"):
        value = config[key]
        if value is not None and (not isinstance(value, str) or not value.strip()):
            raise ValueError(f"{path}: {key} must be a non-empty string")
    if not config["dedup_heap"]:
        raise ValueError(f"{path}: dedup_heap must be a non-empty string")
    for key in NEGATIVE_KEYS:
        if not isinstance(config[key], bool):
            raise ValueError(f"{path}: {key} must be true or false")


def load_config(explicit_path=None):
    path = resolved_config_path(explicit_path)
    config = _defaults()
    if os.path.isfile(path):
        with open(path) as handle:
            data = json.load(handle)
        if not isinstance(data, dict):
            raise ValueError(f"{path}: config must be a JSON object")
        for key in config:
            if key in data:
                config[key] = data[key]
        for positive, negative in POSITIVE_TO_NEGATIVE.items():
            if positive in data:
                config[negative] = not bool(data[positive])
    _validate(config, path)
    return config, path
