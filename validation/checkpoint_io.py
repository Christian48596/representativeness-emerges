"""Safe gzip/JSON checkpoint I/O for image-based structural analyses.

Arrays carry dtype/shape explicitly; no pickle or object execution is used.
Integer-valued mapping keys are restored for support-size lookup.
"""
from pathlib import Path
import gzip
import json
import math
import numpy as np


def _encode(value):
    if isinstance(value, np.ndarray):
        return {'__ndarray__': value.tolist(), 'dtype': value.dtype.str, 'shape': list(value.shape)}
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def _decode(value):
    if set(value)=={'__nonfinite_float__'}:
        token=value['__nonfinite_float__']
        if token not in {'nan','inf','-inf'}:raise ValueError('Invalid nonfinite float token')
        return float(token)
    if '__ndarray__' in value and set(value) == {'__ndarray__', 'dtype', 'shape'}:
        return np.asarray(value['__ndarray__'], dtype=value['dtype']).reshape(value['shape'])
    return {int(k) if k.isdecimal() else k: v for k, v in value.items()}


def save_checkpoint(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Undefined ancillary diagnostics (e.g. constant-field boundary residuals)
    # are preserved without emitting invalid JSON NaN/Infinity literals.
    def prepare(obj):
        if isinstance(obj,np.ndarray):
            return {'__ndarray__':prepare(obj.tolist()),'dtype':obj.dtype.str,'shape':list(obj.shape)}
        if isinstance(obj,np.generic):return prepare(obj.item())
        if isinstance(obj,float) and not math.isfinite(obj):
            return {'__nonfinite_float__':str(obj)}
        if isinstance(obj,dict):return {k:prepare(v) for k,v in obj.items()}
        if isinstance(obj,(list,tuple)):return [prepare(v) for v in obj]
        return obj
    data = json.dumps(prepare(value), separators=(',', ':'), allow_nan=False).encode('utf-8')
    path.write_bytes(gzip.compress(data, mtime=0))


def load_checkpoint(path):
    with gzip.open(path, 'rt', encoding='utf-8') as stream:
        return json.load(stream, object_hook=_decode)


def finite_json(value):
    """Strict JSON: non-finite results become null, interpreted in README."""
    if isinstance(value, np.ndarray):
        return finite_json(value.tolist())
    if isinstance(value, np.generic):
        return finite_json(value.item())
    if isinstance(value, dict):
        return {str(k): finite_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite_json(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_result(path, value):
    Path(path).write_text(json.dumps(finite_json(value), indent=2, allow_nan=False) + '\n', encoding='utf-8')
