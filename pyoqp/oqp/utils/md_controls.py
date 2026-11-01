"""Pure helpers for public ``md(...)`` control provenance and RNG seeds."""

from datetime import date


_PROVENANCE_KEYS = ("common_controls", "common_control_keys")


def explicit_md_control_keys(md_cfg, written_keys=None):
    """Return the ``[md]`` controls the user actually supplied.

    The concise ``md(...)`` call and the Python API record them in
    ``common_control_keys``.  A sectioned ``.inp`` deck has no such marker,
    and without it a plain ``[md] nstep = 2`` was indistinguishable from the
    schema default and was ignored by QM/MM MD (which then ran its own 1000
    steps).  ``written_keys`` are the keys present in the deck as written;
    they are used when there is no marker.
    """
    raw = str(md_cfg.get("common_control_keys", "") or "")
    marked = {key for key in raw.replace(",", " ").split() if key}
    if marked or written_keys is None:
        return marked
    return {str(key) for key in written_keys if str(key) not in _PROVENANCE_KEYS}


def sectioned_md_keys(path):
    """Keys written in the ``[md]`` section of a sectioned deck, or None when
    ``path`` is not such a file (a concise ``.oqp`` input, nothing at all)."""
    import configparser
    import os

    if not path or not os.path.isfile(str(path)):
        return None
    parser = configparser.ConfigParser(interpolation=None, strict=False)
    try:
        with open(str(path), encoding="utf-8") as stream:
            parser.read_file(stream)
    except (configparser.Error, OSError, UnicodeDecodeError):
        return None
    if not parser.has_section("md"):
        return set() if parser.sections() else None
    return set(parser["md"].keys())


def continuation_seed(seed, step):
    """Thermostat seed for a run continued at ``step``.

    A restart builds a new integrator.  Seeding it with the user's seed again
    would replay the random sequence of the first segment in every segment --
    correlated noise, not a continuation.  Mixing in the step gives each
    continuation point its own reproducible stream, in OpenMM's positive
    32-bit range.
    """
    mask = (1 << 64) - 1
    mixed = (int(seed) & mask) ^ (((int(step) + 1) * 0xBF58476D1CE4E5B9) & mask)
    mixed ^= mixed >> 31
    return int(mixed % 2147483646) + 1


def merge_explicit_md_controls(qmmm_cfg, md_cfg):
    """Apply explicit common MD controls without leaking schema defaults."""
    merged = dict(qmmm_cfg)
    explicit = explicit_md_control_keys(md_cfg)
    mapping = (
        ("nstep", "nstep", "n_steps"),
        ("dt", "dt", "timestep"),
        ("temperature", "thermostat_temperature", "temperature"),
        ("thermostat_temperature", "thermostat_temperature", "temperature"),
        ("init_temp", "init_temp", "initial_temperature"),
        ("friction", "thermostat_friction", "friction"),
        ("thermostat_friction", "thermostat_friction", "friction"),
        ("trajectory_file", "trajectory_file", "trajectory_file"),
        ("energy_file", "energy_file", "energy_file"),
        # md(trajectory_interval=N) is the public output cadence; this driver
        # calls it report_interval (frames and log rows share it).  Left
        # unmapped, the request was accepted and every step was written.
        ("trajectory_interval", "trajectory_interval", "report_interval"),
    )
    for provenance_key, md_key, qmmm_key in mapping:
        if provenance_key in explicit and md_key in md_cfg:
            merged[qmmm_key] = md_cfg[md_key]
    if "ensemble" in explicit and "ensemble" in md_cfg:
        # The ensemble is the canonical control and wins.  A legacy thermostat
        # named in an earlier call stays in the provenance list, and letting it
        # decide here would turn a requested npt into nvt -- no barostat, and
        # no message.
        merged["ensemble"] = str(md_cfg["ensemble"]).strip().lower()
    elif "thermostat" in explicit:
        thermostat = str(md_cfg.get("thermostat", "off")).strip().lower()
        merged["ensemble"] = "nvt" if thermostat == "langevin" else "nve"
    return merged


def resolve_output_name(cfg, key, default):
    """An output file name taken from a materialised config section.

    The ``[qmmm]`` output-name defaults are empty strings in the schema, so a
    Runner-materialised config carries the key with an empty value and
    ``cfg.get(key, default)`` hands back ``""`` instead of *default* -- the
    driver then opens ``''`` and dies with ``FileNotFoundError``.  An empty or
    whitespace-only name means unset, exactly as it does on the raw-deck path.
    """
    name = str(cfg.get(key) or "").strip()
    return name or default


def openmm_random_seed(seed, rng_stream, *, today=None):
    """Map the public date/stream seed pair to OpenMM's positive 32-bit seed."""
    base = int(seed)
    if base == 0:
        today = today or date.today()
        base = int(today.strftime("%Y%m%d"))
    stream = int(rng_stream)
    if stream < 0:
        raise ValueError("[md] rng_stream must be non-negative")
    mixed = (
        (base & ((1 << 64) - 1))
        ^ ((stream * 0x9E3779B97F4A7C15) & ((1 << 64) - 1))
    )
    return int(mixed % 2147483646) + 1
