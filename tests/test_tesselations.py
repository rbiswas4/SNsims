import snsims
import pytest


@pytest.mark.skip(reason='not really implemented')
def test_tesselations_abc():

    class NotTile(snsims.Tiling):
        pass
    with pytest.raises(Exception) as e:
            t = NotTile()
    assert e.typename == 'TypeError'
    assert e.value.message == "Can't instantiate abstract class NotTile with abstract methods __init__, area, pointingSequenceForTile, positions, tileIDSequence, tileIDsForSN"
