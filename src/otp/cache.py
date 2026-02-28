import hashlib
import json
import os
import shutil
from pathlib import Path
from typing import Dict, Any, Optional

class CacheManager:
    def __init__(self, runs_dir: Path):
        self.cache_dir = runs_dir / ".cache"
        os.makedirs(self.cache_dir, exist_ok=True)
        # we could also store code_version etc.
        self.version = "v2_0.1.0"
        
    def _generate_hash(self, params: Dict[str, Any]) -> str:
        param_str = json.dumps(params, sort_keys=True)
        combined = f"{self.version}_{param_str}"
        return hashlib.md5(combined.encode()).hexdigest()
        
    def get_cache_path(self, params: Dict[str, Any]) -> str:
        """Return the directory for a specific parameter hash."""
        h = self._generate_hash(params)
        return str(self.cache_dir / h)
        
    def is_cached(self, params: Dict[str, Any]) -> bool:
        """Check if a complete cache exists for these params."""
        cache_path = Path(self.get_cache_path(params))
        manifest_path = cache_path / "manifest.json"
        
        if cache_path.exists() and manifest_path.exists():
            return True
        return False
        
    def load_cache(self, params: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        if not self.is_cached(params):
            return None
            
        cache_path = Path(self.get_cache_path(params))
        manifest_path = cache_path / "manifest.json"
        with open(manifest_path, 'r') as f:
            data = json.load(f)
            
        return {
            "path": str(cache_path),
            "manifest": data
        }
        
    def save_cache(self, params: Dict[str, Any], results: Dict[str, Any], result_files: Dict[str, str]):
        """
        Save results to a new cache directory.
        result_files: { "cas_offinder": "/path/to/cas_output.csv", ... }
        """
        cache_path = Path(self.get_cache_path(params))
        os.makedirs(cache_path, exist_ok=True)
        
        manifest = {
            "version": self.version,
            "params": params,
            "files": {}
        }
        
        # Copy result files into cache directory
        for key, filepath in result_files.items():
            if os.path.exists(filepath):
                dest = cache_path / Path(filepath).name
                shutil.copy2(filepath, dest)
                manifest["files"][key] = dest.name
                
        # Write manifest
        with open(cache_path / "manifest.json", 'w') as f:
            json.dump(manifest, f, indent=2)
