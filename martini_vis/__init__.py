# Find the data directory once.
try:
    from importlib.resources import files, as_file
    import atexit
    from contextlib import ExitStack
except ImportError:
    from pathlib import Path
    DATA_PATH = Path(__file__).parent / 'data'
    del Path
else:
    ref_data = files('martini_vis') / 'data'
    file_manager = ExitStack()
    atexit.register(file_manager.close)
    DATA_PATH = file_manager.enter_context(as_file(ref_data))

    del files, as_file, atexit, ExitStack

from .src.system_reading import system_reading
from .src import index_writer
from .src.molecule_editing import molecule_editor
from .src.topology import topol_writing
