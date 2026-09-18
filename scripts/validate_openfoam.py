"""用独立 VTK OpenFOAM 读取器验证案例；仅为可选测试依赖。"""

import argparse
from pathlib import Path

from vtkmodules.vtkIOGeometry import vtkOpenFOAMReader


def validate(marker):
    errors = []
    reader = vtkOpenFOAMReader()
    reader.AddObserver("ErrorEvent", lambda *_: errors.append("VTK reader error"))
    reader.SetFileName(str(marker.resolve()))
    reader.SetSkipZeroTime(0)
    reader.UpdateInformation()
    for index in range(reader.GetNumberOfCellArrays()):
        reader.SetCellArrayStatus(reader.GetCellArrayName(index), 1)
    times = reader.GetTimeValues()
    assert times is not None and times.GetNumberOfValues() > 0, marker
    assert times.GetValue(0) == 0, (marker, "Missing initial time")
    count = None
    for index in range(times.GetNumberOfValues()):
        time = times.GetValue(index)
        reader.UpdateTimeStep(time)
        assert not errors, (marker, errors)
        blocks = reader.GetOutput()
        mesh = blocks.GetBlock(0)
        assert mesh is not None and mesh.IsA("vtkUnstructuredGrid"), marker
        assert mesh.GetNumberOfCells() > 0, marker
        if count is None:
            count = mesh.GetNumberOfCells()
        assert mesh.GetNumberOfCells() == count, marker
        data = mesh.GetCellData()
        for name in ("U", "rho", "p", "porosity"):
            array = data.GetArray(name)
            assert array is not None and array.GetNumberOfTuples() == count, (marker, name)
        assert data.GetArray("U").GetNumberOfComponents() == 3, marker
        a, b = data.GetArray("alpha.A"), data.GetArray("alpha.B")
        if a is not None:
            assert b is not None and data.GetArray("pBulk") is not None, marker
            for cell in range(count):
                assert abs(a.GetValue(cell) + b.GetValue(cell) - 1) < 1e-6, marker
    print(f"PASS {marker}: {count} cells, {times.GetNumberOfValues()} times")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path)
    args = parser.parse_args()
    files = []
    for path in args.paths:
        files.extend([path] if path.is_file() else sorted(path.rglob("*.foam")))
    if not files:
        raise SystemExit("No .foam files found")
    for path in files:
        validate(path)
