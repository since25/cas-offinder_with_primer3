import os
from pathlib import Path
from otp.cache import CacheManager

def test_cache_hashing(tmp_path):
    mgr = CacheManager(tmp_path)
    
    params1 = {"spacer": "AAA", "pam": "NGG", "mismatches": 0}
    params2 = {"spacer": "AAA", "mismatches": 0, "pam": "NGG"}
    # The hash should be identical regardless of key order
    h1 = mgr._generate_hash(params1)
    h2 = mgr._generate_hash(params2)
    
    assert h1 == h2

def test_cache_save_load(tmp_path):
    mgr = CacheManager(tmp_path)
    params = {"query": "test"}
    
    assert mgr.is_cached(params) is False
    
    mgr.save_cache(params, [{"a": 1}], {})
    
    assert mgr.is_cached(params) is True
    data = mgr.load_cache(params)
    assert data is not None
    assert "version" in data["manifest"]
