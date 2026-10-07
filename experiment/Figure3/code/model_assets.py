"""Resolve model assets locally or download verified public checkpoints to cache."""
import hashlib
from pathlib import Path
import urllib.request

from huggingface_hub import hf_hub_download


HF_REPO_ID = "aer0vane/reproduce_deepsea"
OURS_CHECKPOINT_NAME = "best_model_FINAL_EPOCH53.pth"
OURS_CHECKPOINT_HF_PATH = f"training_checkpoints/{OURS_CHECKPOINT_NAME}"
PRETRAINED_PREDICT_NAME = "deepsea_predict.pth"
PRETRAINED_VARIANT_NAME = "deepsea_variant_effects.pth"
PRETRAINED_PREDICT_URL = (
    "https://zenodo.org/records/1466993/files/"
    "deepsea_predict.pth?download=1"
)
PRETRAINED_VARIANT_URL = (
    "https://zenodo.org/records/1466993/files/"
    "deepsea_variant_effects.pth?download=1"
)
PRETRAINED_MD5 = {
    PRETRAINED_PREDICT_NAME: "89e640bf6bdbe1ff165f484d9796efc7",
    PRETRAINED_VARIANT_NAME: "35956ab9c28960b5a3693f470fe980c1",
}


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _download_verified_checkpoint(name: str, url: str) -> Path:
    cache_dir = Path.home() / ".cache" / "deepsea-fuxian3"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached = cache_dir / name
    if cached.is_file() and _md5(cached) == PRETRAINED_MD5[name]:
        return cached
    partial = cached.with_suffix(cached.suffix + ".part")
    try:
        urllib.request.urlretrieve(url, partial)
        actual = _md5(partial)
        expected = PRETRAINED_MD5[name]
        if actual != expected:
            raise ValueError(f"DeepSEA {name} MD5 mismatch: {actual} != {expected}")
        partial.replace(cached)
    finally:
        partial.unlink(missing_ok=True)
    return cached


def resolve_pretrained_predict_checkpoint(project_root: str | Path) -> Path:
    """Resolve the released sequence-prediction checkpoint (919 chromatin outputs)."""
    local_path = Path(project_root) / "models" / PRETRAINED_PREDICT_NAME
    if local_path.is_file():
        return local_path
    try:
        downloaded = Path(hf_hub_download(
            repo_id=HF_REPO_ID,
            filename=f"models/{PRETRAINED_PREDICT_NAME}",
            repo_type="model",
        ))
        if _md5(downloaded) == PRETRAINED_MD5[PRETRAINED_PREDICT_NAME]:
            return downloaded
    except Exception:
        pass
    # Zenodo is the authoritative upstream source and public fallback.
    return _download_verified_checkpoint(PRETRAINED_PREDICT_NAME, PRETRAINED_PREDICT_URL)


def resolve_pretrained_variant_effects_checkpoint(project_root: str | Path) -> Path:
    """Resolve the released variant-effects checkpoint (not the sequence-only model)."""
    local_path = Path(project_root) / "models" / PRETRAINED_VARIANT_NAME
    if local_path.is_file():
        return local_path
    return _download_verified_checkpoint(PRETRAINED_VARIANT_NAME, PRETRAINED_VARIANT_URL)


def resolve_ours_checkpoint(project_root: str | Path) -> Path:
    """Return a local trained checkpoint, downloading the published one if absent."""
    local_path = Path(project_root) / "training" / "checkpoints" / OURS_CHECKPOINT_NAME
    if local_path.is_file():
        return local_path
    cached_path = hf_hub_download(
        repo_id=HF_REPO_ID,
        filename=OURS_CHECKPOINT_HF_PATH,
        repo_type="model",
    )
    return Path(cached_path)
