import itsxpress.main
import itsxpress.definitions
import itsxpress.ITSposition
import itsxpress.Dedup
import itsxpress.SeqSample
import itsxpress.SeqSampleNotPaired
import itsxpress.SeqSamplePaired

try:
    import itsxpress.q2_itsxpress
    import itsxpress.plugin_setup
except ModuleNotFoundError as e:
    #logging
    print("{}.Could not initialize the Qiime plugin portion of ITSxpress. Command line ITSxpress will still work normally. If you wish to use the Qiime2 ITSxpress plugin, you need to install Qiime2 first into your environment.\n".format(e))
    pass

__all__ = ["main", "definitions"]

from . import _version
__version__ = _version.get_versions()['version']
