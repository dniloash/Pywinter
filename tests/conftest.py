from pathlib import Path
import pytest


@pytest.fixture
def ex_data_1(request: pytest.FixtureRequest):
    return (
        Path(request.config.rootpath)
        / "examples/data_example_01/MERRA2_400.tavg1_2d_slv_Nx.20150105.SUB.nc"
    )
