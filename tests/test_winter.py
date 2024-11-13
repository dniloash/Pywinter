from pathlib import Path
import pytest

import numpy as np
from netCDF4 import Dataset, num2date
import pywinter.winter as pyw


def test_geo0():
    stlat, stlon = 0.0, 0.0
    dlat, dlon = 0.5, 0.5
    geo = pyw.Geo0(stlat, stlon, dlat, dlon)
    assert geo is not None


V2D = (
    "PSFC",
    "PMSL",
    "SKINTEMP",
    "SOILHGT",
    "TT",
    "RH",
    "SPECHUMD",
    "UU",
    "VV",
    "LANDSEA",
    "SST",
    "SEAICE",
    "SNOW",
    "TAVGSFC",
)


def test_V2d(ex_data_1):
    netfile = Dataset(ex_data_1, "r")
    var = netfile.variables["T2M"][0, :, :]
    var3d = netfile.variables["T2M"][:, :, :]

    w_var = pyw.V2d(V2D[0], var)
    assert w_var.name == V2D[0]
    assert w_var.des == ""
    assert w_var.uni == ""
    assert w_var.lev == ""

    des, uni, lev = w_var.idvar()

    w_var = pyw.V2d(V2D[0], var, des, uni, lev)
    assert w_var.name == V2D[0]
    assert w_var.des == des
    assert w_var.uni == uni
    assert w_var.lev == lev

    w_var = pyw.V2d("RANDOMVARNAME", var, des, uni, lev)
    with pytest.raises(Exception):
        w_var.idvar()

    with pytest.raises(Exception):
        pyw.V2d(V2D[0], var3d)

    netfile.close()


def test_io(tmp_path, ex_data_1: Path):
    netfile = Dataset(ex_data_1, "r")

    temp2m = netfile.variables["T2M"][:, :, :]
    lat = netfile.variables["lat"][:]
    lon = netfile.variables["lon"][:]

    dlat = np.abs(lat[1] - lat[0])
    dlon = np.abs(lon[1] - lon[0])

    time = netfile.variables["time"]
    timed1 = num2date(time[:], units=time.units, calendar=time.calendar)
    timed2 = [str(i)[:13].replace(" ", "_") for i in timed1]

    geo = pyw.Geo0(lat[0], lon[0], dlat, dlon)

    for i in range(len(timed2)):
        winter_t2 = pyw.V2d("TT", temp2m[i, :, :])
        total_var = [winter_t2]
        pyw.cinter("FILETEMP", timed2[i], geo, total_var, str(tmp_path))
    inter_files = list(tmp_path.glob("FILETEMP*"))
    assert len(inter_files) == len(timed2)

    for f in inter_files:
        ff = pyw.rinter(f)
        assert len(ff.keys()) >= 0

    netfile.close()
