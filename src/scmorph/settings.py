from tempfile import TemporaryDirectory

import scanpy

settings = scanpy.settings
settings.cachedir = TemporaryDirectory()
