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
except Exception as e:
    print("could not initialize the Qiime plugin portion of ITSxpress. Command line ITSxpress should work normally")
    print(e)
    pass

__all__ = ["main", "definitions"]
