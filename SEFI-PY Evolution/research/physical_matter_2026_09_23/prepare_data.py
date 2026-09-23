"""Unpack research data locally; never overwrite different files."""
from pathlib import Path
import zipfile
root=Path(__file__).resolve().parent
with zipfile.ZipFile(root/'reference_profiles.zip') as archive:
    for item in archive.infolist():
        target=(root/item.filename).resolve()
        if not target.is_relative_to(root) or not item.filename.startswith('work/'):
            raise ValueError('Unexpected archive path')
        data=archive.read(item)
        if target.exists() and target.read_bytes()!=data:
            raise FileExistsError('Refusing to overwrite changed research data: '+str(target))
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(data)
print('Reference data ready in this isolated research folder.')
