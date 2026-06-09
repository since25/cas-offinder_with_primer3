import gzip

import pytest

from scripts.download_genomes import decompress_gzip


def test_decompress_gzip_does_not_leave_final_file_on_bad_archive(tmp_path):
    gz_path = tmp_path / "bad.fa.gz"
    out_path = tmp_path / "genome.fa"

    with gzip.open(gz_path, "wb") as handle:
        handle.write(b"ACGT" * 1024)

    data = gz_path.read_bytes()
    gz_path.write_bytes(data[:-8] + b"BADCRC!!")

    with pytest.raises(gzip.BadGzipFile):
        decompress_gzip(gz_path, out_path)

    assert not out_path.exists()
    assert not out_path.with_suffix(out_path.suffix + ".tmp").exists()
